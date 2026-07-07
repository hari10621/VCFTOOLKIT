"""VCF validation helpers."""

from __future__ import annotations

from pathlib import Path

from .exceptions import VCFValidationError
from .utils import open_text_stream

REQUIRED_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


class VCFValidator:
    """Validates VCF files without loading them into memory."""

    def validate_path(self, path: str | Path) -> None:
        """Validate file existence, extension, header, and first data rows."""

        path = Path(path)
        if not path.exists():
            raise VCFValidationError(f"VCF file does not exist: {path}")
        if path.suffix.lower() not in {".vcf", ".gz"} and not path.name.lower().endswith(".vcf.gz"):
            raise VCFValidationError("Unsupported file type. Expected .vcf or .vcf.gz")
        if path.stat().st_size == 0:
            raise VCFValidationError("VCF file is empty")

        found_fileformat = False
        found_columns = False
        data_checked = 0
        try:
            with open_text_stream(path) as handle:
                for line_number, line in enumerate(handle, start=1):
                    line = line.rstrip("\n")
                    if line.startswith("##fileformat=VCF"):
                        found_fileformat = True
                        continue
                    if line.startswith("#CHROM"):
                        columns = line.split("\t")
                        if columns[:8] != REQUIRED_COLUMNS:
                            raise VCFValidationError(
                                f"Invalid VCF column header at line {line_number}: expected first 8 columns"
                            )
                        found_columns = True
                        continue
                    if line.startswith("#"):
                        continue
                    if found_columns and line:
                        columns = line.split("\t")
                        if len(columns) < 8:
                            raise VCFValidationError(
                                f"Malformed variant line {line_number}: expected at least 8 tab-separated columns"
                            )
                        self._validate_variant_columns(columns, line_number)
                        data_checked += 1
                        if data_checked >= 10:
                            break
        except UnicodeError as exc:
            raise VCFValidationError(f"Unable to decode VCF file: {exc}") from exc
        except OSError as exc:
            raise VCFValidationError(f"Unable to read VCF file: {exc}") from exc

        if not found_fileformat:
            raise VCFValidationError("Missing required ##fileformat=VCF header")
        if not found_columns:
            raise VCFValidationError("Missing required #CHROM column header")

    def _validate_variant_columns(self, columns: list[str], line_number: int) -> None:
        try:
            position = int(columns[1])
        except ValueError as exc:
            raise VCFValidationError(f"Invalid POS value at line {line_number}: {columns[1]}") from exc
        if position < 1:
            raise VCFValidationError(f"Invalid POS value at line {line_number}: must be positive")
        if not columns[3] or columns[3] == ".":
            raise VCFValidationError(f"Missing REF value at line {line_number}")
        if not columns[4] or columns[4] == ".":
            raise VCFValidationError(f"Missing ALT value at line {line_number}")
