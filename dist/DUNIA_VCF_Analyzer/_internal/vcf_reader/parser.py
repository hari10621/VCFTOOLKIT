"""High-performance streaming VCF reader."""

from __future__ import annotations

import logging
import threading
import time
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

from .exceptions import VCFParseError, VCFReaderStateError, VCFValidationError
from .exporter import VCFExporter
from .filters import filter_stream, search_variant
from .models import ParserBackend, ProgressInfo, VCFHeader, VCFVariant, VariantFilter, VariantStatistics
from .statistics import StatisticsCollector
from .utils import (
    META_PATTERN,
    STRUCTURED_META_PATTERN,
    detect_backend,
    get_file_size,
    open_text_stream,
    parse_filter,
    parse_info,
    parse_key_value_body,
    parse_samples,
)
from .validator import VCFValidator

logger = logging.getLogger(__name__)


class VCFReader:
    """Streaming reader for large VCF and VCF.GZ files."""

    def __init__(self, path: str | Path | None = None, backend: ParserBackend | None = None) -> None:
        self.path: Path | None = None
        self.backend = detect_backend(backend)
        self.header: VCFHeader | None = None
        self.statistics = StatisticsCollector()
        self.validator = VCFValidator()
        self.exporter = VCFExporter()
        self._file_size = 0
        self._file_position = 0
        self._variant_count = 0
        self._started_at: float | None = None
        self._cancel_event = threading.Event()
        self._pause_event = threading.Event()
        self._pause_event.set()
        self._worker_thread: threading.Thread | None = None
        self._lock = threading.Lock()
        if path:
            self.open_file(path)

    def open_file(self, path: str | Path) -> None:
        """Open a VCF path and select the best available backend."""

        self.path = Path(path)
        self._file_size = get_file_size(self.path)
        self.header = None
        self.statistics = StatisticsCollector()
        self._variant_count = 0
        self._file_position = 0
        self._cancel_event.clear()
        self._pause_event.set()
        logger.info("File opened: %s using backend=%s", self.path, self.backend.value)

    def validate(self) -> bool:
        """Validate the currently opened file."""

        self._require_path()
        self.validator.validate_path(self.path)
        logger.info("VCF validation passed: %s", self.path)
        return True

    def read_header(self) -> VCFHeader:
        """Read and parse VCF header metadata."""

        self._require_path()
        header = VCFHeader()
        try:
            with open_text_stream(self.path) as handle:
                for line in handle:
                    stripped = line.rstrip("\n")
                    if stripped.startswith("##"):
                        header.raw_lines.append(stripped)
                        self._parse_header_line(stripped, header)
                    elif stripped.startswith("#CHROM"):
                        header.column_line = stripped
                        columns = stripped.split("\t")
                        if len(columns) < 8:
                            raise VCFValidationError("Invalid #CHROM header: expected at least 8 columns")
                        header.samples = columns[9:] if len(columns) > 9 else []
                        break
                    elif stripped:
                        raise VCFValidationError("VCF data appeared before #CHROM header")
        except OSError as exc:
            raise VCFValidationError(f"Unable to read header: {exc}") from exc
        if not header.version:
            raise VCFValidationError("Missing ##fileformat VCF header")
        if not header.column_line:
            raise VCFValidationError("Missing #CHROM VCF column header")
        self.header = header
        logger.info("Header parsed: version=%s samples=%d", header.version, len(header.samples))
        return header

    def stream_variants(self, update_statistics: bool = True) -> Iterator[VCFVariant]:
        """Yield variants one by one without retaining the whole file."""

        self._require_path()
        if self.header is None:
            self.read_header()
        self._started_at = self._started_at or time.monotonic()
        if self.backend == ParserBackend.CYVCF2:
            yield from self._stream_cyvcf2(update_statistics)
        elif self.backend == ParserBackend.PYSAM:
            yield from self._stream_pysam(update_statistics)
        else:
            yield from self._stream_manual(update_statistics)

    def filter_variants(self, criteria: VariantFilter | None = None, **kwargs: Any) -> Iterator[VCFVariant]:
        """Yield variants matching runtime criteria."""

        criteria = criteria or VariantFilter(**kwargs)
        yield from filter_stream(self.stream_variants(), criteria)

    def search(self, limit: int | None = None, **query: Any) -> Iterator[VCFVariant]:
        """Search variants by gene, chromosome, position, ID, reference, or alternate."""

        matches = 0
        for variant in self.stream_variants(update_statistics=False):
            if search_variant(variant, **query):
                yield variant
                matches += 1
                if limit is not None and matches >= limit:
                    return

    def get_statistics(self, refresh: bool = False) -> VariantStatistics:
        """Return streaming statistics, optionally rescanning the file."""

        if refresh or self.statistics.stats.total_variants == 0:
            self.statistics = StatisticsCollector()
            for _ in self.stream_variants(update_statistics=True):
                if self._cancel_event.is_set():
                    break
        return self.statistics.snapshot()

    def export(
        self,
        output_path: str | Path,
        export_format: str | None = None,
        criteria: VariantFilter | None = None,
    ) -> int:
        """Export all or filtered variants to CSV, TSV, JSON, Excel, or VCF."""

        output_path = Path(output_path)
        export_format = export_format or output_path.suffix.lstrip(".")
        variants = self.filter_variants(criteria) if criteria else self.stream_variants(update_statistics=False)
        written = self.exporter.export(variants, output_path, export_format, self.header)
        logger.info("Export completed: %s variants=%d", output_path, written)
        return written

    def run_in_background(
        self,
        on_variant: Callable[[VCFVariant], None] | None = None,
        on_complete: Callable[[VariantStatistics], None] | None = None,
        on_error: Callable[[Exception], None] | None = None,
    ) -> threading.Thread:
        """Start parsing in a background worker thread."""

        if self._worker_thread and self._worker_thread.is_alive():
            raise VCFReaderStateError("Background parser is already running")
        self._cancel_event.clear()
        self._pause_event.set()

        def worker() -> None:
            try:
                for variant in self.stream_variants(update_statistics=True):
                    if on_variant:
                        on_variant(variant)
                    if self._cancel_event.is_set():
                        logger.info("Parsing cancelled")
                        return
                if on_complete:
                    on_complete(self.statistics.snapshot())
            except Exception as exc:  # pragma: no cover - callback plumbing
                logger.exception("Background parsing failed")
                if on_error:
                    on_error(exc)

        self._worker_thread = threading.Thread(target=worker, name="VCFReaderWorker", daemon=True)
        self._worker_thread.start()
        return self._worker_thread

    def pause(self) -> None:
        """Pause background or foreground streaming at the next record."""

        self._pause_event.clear()
        logger.info("Parsing paused")

    def resume(self) -> None:
        """Resume paused parsing."""

        self._pause_event.set()
        logger.info("Parsing resumed")

    def cancel(self) -> None:
        """Request cancellation."""

        self._cancel_event.set()
        self._pause_event.set()
        logger.info("Parsing cancellation requested")

    def progress(self) -> ProgressInfo:
        """Return current parsing progress."""

        elapsed = time.monotonic() - self._started_at if self._started_at else 0.0
        speed = self._variant_count / elapsed if elapsed > 0 else 0.0
        percent = (self._file_position / self._file_size * 100) if self._file_size else 0.0
        remaining = None
        if speed > 0 and percent > 0:
            remaining = elapsed * ((100 - percent) / percent)
        return ProgressInfo(
            file_position=self._file_position,
            file_size=self._file_size,
            percent_complete=min(percent, 100.0),
            estimated_time_remaining=remaining,
            current_variant_count=self._variant_count,
            variants_per_second=speed,
            elapsed_time=elapsed,
        )

    def close(self) -> None:
        """Cancel active work and release reader state."""

        self.cancel()
        if self._worker_thread and self._worker_thread.is_alive():
            self._worker_thread.join(timeout=5)
        self._worker_thread = None
        logger.info("VCF reader closed")

    def _stream_manual(self, update_statistics: bool) -> Iterator[VCFVariant]:
        samples = self.header.samples if self.header else []
        try:
            with open_text_stream(self.path) as handle:
                for line_number, line in enumerate(handle, start=1):
                    self._file_position = self._safe_tell(handle)
                    if line.startswith("#"):
                        continue
                    if self._cancel_event.is_set():
                        return
                    self._pause_event.wait()
                    stripped = line.rstrip("\n")
                    if not stripped:
                        continue
                    try:
                        variant = self._parse_variant_line(stripped, samples)
                    except VCFParseError:
                        logger.warning("Malformed variant skipped at line %d", line_number)
                        raise
                    self._record_variant(variant, update_statistics)
                    yield variant
        except OSError as exc:
            raise VCFParseError(f"Unable to stream variants: {exc}") from exc

    def _stream_cyvcf2(self, update_statistics: bool) -> Iterator[VCFVariant]:
        try:
            from cyvcf2 import VCF
        except ImportError:
            self.backend = ParserBackend.MANUAL
            yield from self._stream_manual(update_statistics)
            return

        for record in VCF(str(self.path)):
            if self._cancel_event.is_set():
                return
            self._pause_event.wait()
            variant = VCFVariant(
                chromosome=record.CHROM,
                position=int(record.POS),
                id=None if record.ID in {None, "."} else record.ID,
                reference=record.REF,
                alternate=list(record.ALT or []),
                qual=float(record.QUAL) if record.QUAL is not None else None,
                filter=["PASS"] if record.FILTER in {None, "PASS"} else str(record.FILTER).split(";"),
                info=dict(record.INFO),
                format_keys=[],
                sample_values={},
            )
            self._record_variant(variant, update_statistics)
            yield variant

    def _stream_pysam(self, update_statistics: bool) -> Iterator[VCFVariant]:
        try:
            import pysam
        except ImportError:
            self.backend = ParserBackend.MANUAL
            yield from self._stream_manual(update_statistics)
            return

        with pysam.VariantFile(str(self.path)) as handle:
            for record in handle.fetch():
                if self._cancel_event.is_set():
                    return
                self._pause_event.wait()
                variant = VCFVariant(
                    chromosome=record.chrom,
                    position=int(record.pos),
                    id=None if record.id in {None, "."} else record.id,
                    reference=record.ref,
                    alternate=list(record.alts or []),
                    qual=float(record.qual) if record.qual is not None else None,
                    filter=list(record.filter.keys()) or ["PASS"],
                    info=dict(record.info),
                    format_keys=list(record.format.keys()),
                    sample_values={sample: dict(values) for sample, values in record.samples.items()},
                )
                self._record_variant(variant, update_statistics)
                yield variant

    def _parse_header_line(self, line: str, header: VCFHeader) -> None:
        if line.startswith("##fileformat="):
            header.version = line.split("=", 1)[1]
            header.metadata.setdefault("fileformat", []).append(header.version)
            return
        structured = STRUCTURED_META_PATTERN.match(line)
        if structured:
            key = structured.group("key")
            values = parse_key_value_body(structured.group("body"))
            identifier = values.get("ID")
            if identifier and key == "INFO":
                header.info_definitions[identifier] = values
            elif identifier and key == "FILTER":
                header.filter_definitions[identifier] = values
            elif identifier and key == "FORMAT":
                header.format_definitions[identifier] = values
            return
        match = META_PATTERN.match(line)
        if not match:
            return
        key = match.group("key")
        value = match.group("value")
        header.metadata.setdefault(key, []).append(value)
        if key.lower() == "reference":
            header.reference = value
        elif key.lower() == "contig":
            parsed = parse_key_value_body(value.strip("<>"))
            if "ID" in parsed:
                header.contigs.append(parsed["ID"])

    def _parse_variant_line(self, line: str, samples: list[str]) -> VCFVariant:
        columns = line.split("\t")
        if len(columns) < 8:
            raise VCFParseError("Malformed variant: expected at least 8 columns")
        try:
            position = int(columns[1])
        except ValueError as exc:
            raise VCFParseError(f"Invalid POS value: {columns[1]}") from exc
        qual = None if columns[5] == "." else float(columns[5])
        format_keys, sample_values = parse_samples(columns[8] if len(columns) > 8 else "", columns[9:], samples)
        return VCFVariant(
            chromosome=columns[0],
            position=position,
            id=None if columns[2] == "." else columns[2],
            reference=columns[3],
            alternate=[] if columns[4] == "." else columns[4].split(","),
            qual=qual,
            filter=parse_filter(columns[6]),
            info=parse_info(columns[7]),
            format_keys=format_keys,
            sample_values=sample_values,
            raw_line=line,
        )

    def _record_variant(self, variant: VCFVariant, update_statistics: bool) -> None:
        with self._lock:
            self._variant_count += 1
        if update_statistics:
            self.statistics.update(variant)
        if self._variant_count % 100000 == 0:
            logger.info("Variants processed: %d", self._variant_count)

    def _safe_tell(self, handle: Any) -> int:
        try:
            raw = getattr(handle, "buffer", None)
            if raw is not None:
                return int(raw.tell())
            return int(handle.tell())
        except (OSError, ValueError):
            return self._file_position

    def _require_path(self) -> None:
        if self.path is None:
            raise VCFReaderStateError("No VCF file has been opened")
