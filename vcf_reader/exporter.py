"""Streaming exporters for VCF variants."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from .exceptions import VCFExportError
from .models import VCFHeader, VCFVariant


class VCFExporter:
    """Exports variant streams to common data formats."""

    FIELDNAMES = [
        "chromosome",
        "position",
        "id",
        "reference",
        "alternate",
        "qual",
        "filter",
        "info",
        "format",
        "samples",
    ]

    def export(
        self,
        variants: Iterable[VCFVariant],
        output_path: str | Path,
        export_format: str,
        header: VCFHeader | None = None,
    ) -> int:
        """Export a variant stream and return the number of written variants."""

        output_path = Path(output_path)
        export_format = export_format.lower().lstrip(".")
        try:
            if export_format == "csv":
                return self._delimited(variants, output_path, ",")
            if export_format == "tsv":
                return self._delimited(variants, output_path, "\t")
            if export_format == "json":
                return self._json(variants, output_path)
            if export_format in {"xlsx", "excel"}:
                return self._excel(variants, output_path)
            if export_format == "vcf":
                return self._vcf(variants, output_path, header)
        except OSError as exc:
            raise VCFExportError(f"Export failed: {exc}") from exc
        raise VCFExportError(f"Unsupported export format: {export_format}")

    def _delimited(self, variants: Iterable[VCFVariant], output_path: Path, delimiter: str) -> int:
        count = 0
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.FIELDNAMES, delimiter=delimiter)
            writer.writeheader()
            for variant in variants:
                writer.writerow(self._row(variant))
                count += 1
        return count

    def _json(self, variants: Iterable[VCFVariant], output_path: Path) -> int:
        count = 0
        with output_path.open("w", encoding="utf-8") as handle:
            handle.write("[\n")
            first = True
            for variant in variants:
                if not first:
                    handle.write(",\n")
                json.dump(self._row(variant), handle, ensure_ascii=False)
                first = False
                count += 1
            handle.write("\n]\n")
        return count

    def _excel(self, variants: Iterable[VCFVariant], output_path: Path) -> int:
        try:
            from openpyxl import Workbook
        except ImportError as exc:
            raise VCFExportError("Excel export requires openpyxl to be installed") from exc

        workbook = Workbook(write_only=True)
        sheet = workbook.create_sheet("variants")
        sheet.append(self.FIELDNAMES)
        count = 0
        for variant in variants:
            row = self._row(variant)
            sheet.append([row[name] for name in self.FIELDNAMES])
            count += 1
        workbook.save(output_path)
        return count

    def _vcf(self, variants: Iterable[VCFVariant], output_path: Path, header: VCFHeader | None) -> int:
        count = 0
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            if header:
                for line in header.raw_lines:
                    handle.write(line.rstrip("\n") + "\n")
                if header.column_line:
                    handle.write(header.column_line.rstrip("\n") + "\n")
            else:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for variant in variants:
                handle.write(self._vcf_line(variant) + "\n")
                count += 1
        return count

    def _row(self, variant: VCFVariant) -> dict[str, Any]:
        return {
            "chromosome": variant.chromosome,
            "position": variant.position,
            "id": variant.id,
            "reference": variant.reference,
            "alternate": ",".join(variant.alternate),
            "qual": variant.qual,
            "filter": ";".join(variant.filter) if variant.filter else ".",
            "info": json.dumps(variant.info, ensure_ascii=False, separators=(",", ":")),
            "format": ":".join(variant.format_keys),
            "samples": json.dumps(variant.sample_values, ensure_ascii=False, separators=(",", ":")),
        }

    def _vcf_line(self, variant: VCFVariant) -> str:
        info = self._format_info(variant.info)
        filt = ";".join(variant.filter) if variant.filter else "."
        qual = "." if variant.qual is None else f"{variant.qual:g}"
        fields = [
            variant.chromosome,
            str(variant.position),
            variant.id or ".",
            variant.reference,
            ",".join(variant.alternate),
            qual,
            filt,
            info,
        ]
        if variant.format_keys:
            fields.append(":".join(variant.format_keys))
            for sample in variant.sample_values.values():
                fields.append(":".join("." if sample.get(key) is None else str(sample.get(key)) for key in variant.format_keys))
        return "\t".join(fields)

    def _format_info(self, info: dict[str, Any]) -> str:
        if not info:
            return "."
        parts = []
        for key, value in info.items():
            if value is True:
                parts.append(key)
            elif isinstance(value, list):
                parts.append(f"{key}={','.join(str(item) for item in value)}")
            else:
                parts.append(f"{key}={value}")
        return ";".join(parts)
