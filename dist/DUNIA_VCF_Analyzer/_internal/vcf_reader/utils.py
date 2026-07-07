"""Utility helpers for streaming VCF processing."""

from __future__ import annotations

import gzip
import importlib.util
import os
import re
from pathlib import Path
from typing import Any, TextIO

from .models import ParserBackend

META_PATTERN = re.compile(r"##(?P<key>[^=]+)=(?P<value>.*)")
STRUCTURED_META_PATTERN = re.compile(r"##(?P<key>INFO|FILTER|FORMAT)=<(?P<body>.*)>")


def detect_backend(preferred: ParserBackend | None = None) -> ParserBackend:
    """Return the best available parser backend."""

    if preferred:
        return preferred
    if importlib.util.find_spec("cyvcf2"):
        return ParserBackend.CYVCF2
    if importlib.util.find_spec("pysam"):
        return ParserBackend.PYSAM
    return ParserBackend.MANUAL


def open_text_stream(path: str | Path) -> TextIO:
    """Open a plain or gzip-compressed VCF as a buffered text stream."""

    path = Path(path)
    try:
        with path.open("rb") as probe:
            magic = probe.read(2)
        if magic == b"\x1f\x8b" or path.suffix.lower() == ".gz":
            return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
        return path.open("r", encoding="utf-8", errors="replace")
    except OSError:
        raise


def get_file_size(path: str | Path) -> int:
    """Return file size in bytes."""

    return os.path.getsize(path)


def parse_key_value_body(body: str) -> dict[str, str]:
    """Parse comma-separated VCF angle-bracket metadata into a dictionary."""

    result: dict[str, str] = {}
    token = []
    in_quote = False
    parts: list[str] = []
    for char in body:
        if char == '"':
            in_quote = not in_quote
        if char == "," and not in_quote:
            parts.append("".join(token))
            token = []
        else:
            token.append(char)
    if token:
        parts.append("".join(token))

    for part in parts:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        result[key.strip()] = value.strip().strip('"')
    return result


def parse_info(value: str) -> dict[str, Any]:
    """Parse a VCF INFO field."""

    if value in {"", "."}:
        return {}
    info: dict[str, Any] = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" not in item:
            info[item] = True
            continue
        key, raw = item.split("=", 1)
        values = raw.split(",")
        parsed = [_coerce_scalar(v) for v in values]
        info[key] = parsed[0] if len(parsed) == 1 else parsed
    return info


def parse_filter(value: str) -> list[str]:
    """Parse the FILTER column."""

    if value in {"", ".", "PASS"}:
        return [] if value in {"", "."} else ["PASS"]
    return value.split(";")


def parse_samples(format_value: str, sample_values: list[str], samples: list[str]) -> tuple[list[str], dict[str, dict[str, Any]]]:
    """Parse FORMAT and sample columns."""

    if not format_value or format_value == ".":
        return [], {}
    keys = format_value.split(":")
    parsed: dict[str, dict[str, Any]] = {}
    for sample, raw_values in zip(samples, sample_values, strict=False):
        values = raw_values.split(":")
        parsed[sample] = {
            key: _coerce_scalar(values[index]) if index < len(values) else None
            for index, key in enumerate(keys)
        }
    return keys, parsed


def _coerce_scalar(value: str) -> Any:
    if value == ".":
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def first_numeric(values: list[Any] | tuple[Any, ...] | Any) -> float | None:
    """Return the first numeric value from a scalar or list-like value."""

    if isinstance(values, (list, tuple)):
        for value in values:
            found = first_numeric(value)
            if found is not None:
                return found
        return None
    if isinstance(values, (int, float)):
        return float(values)
    return None
