"""Typed models used by the VCF reader."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class ParserBackend(str, Enum):
    """Parser implementation selected at runtime."""

    CYVCF2 = "cyvcf2"
    PYSAM = "pysam"
    MANUAL = "manual"


@dataclass(slots=True)
class VCFHeader:
    """Parsed VCF header metadata."""

    version: str | None = None
    reference: str | None = None
    info_definitions: dict[str, dict[str, str]] = field(default_factory=dict)
    filter_definitions: dict[str, dict[str, str]] = field(default_factory=dict)
    format_definitions: dict[str, dict[str, str]] = field(default_factory=dict)
    samples: list[str] = field(default_factory=list)
    metadata: dict[str, list[str]] = field(default_factory=dict)
    contigs: list[str] = field(default_factory=list)
    raw_lines: list[str] = field(default_factory=list)
    column_line: str | None = None


@dataclass(slots=True)
class VCFVariant:
    """A normalized variant record."""

    chromosome: str
    position: int
    id: str | None
    reference: str
    alternate: list[str]
    qual: float | None
    filter: list[str]
    info: dict[str, Any]
    format_keys: list[str] = field(default_factory=list)
    sample_values: dict[str, dict[str, Any]] = field(default_factory=dict)
    raw_line: str | None = None

    @property
    def variant_id(self) -> str | None:
        """Alias for callers that prefer explicit naming."""

        return self.id


@dataclass(slots=True)
class VariantFilter:
    """Runtime filtering options for streaming variants."""

    chromosome: str | None = None
    gene: str | None = None
    position_start: int | None = None
    position_end: int | None = None
    min_qual: float | None = None
    min_depth: int | None = None
    min_allele_frequency: float | None = None
    max_allele_frequency: float | None = None
    variant_type: str | None = None
    pass_only: bool = False
    info_fields: dict[str, Any] = field(default_factory=dict)
    sample: str | None = None


@dataclass(slots=True)
class ProgressInfo:
    """Point-in-time progress data for a running parser."""

    file_position: int
    file_size: int
    percent_complete: float
    estimated_time_remaining: float | None
    current_variant_count: int
    variants_per_second: float
    elapsed_time: float


@dataclass(slots=True)
class VariantStatistics:
    """Streaming aggregate statistics."""

    total_variants: int = 0
    snp_count: int = 0
    insertion_count: int = 0
    deletion_count: int = 0
    structural_variant_count: int = 0
    transition_count: int = 0
    transversion_count: int = 0
    chromosome_distribution: dict[str, int] = field(default_factory=dict)
    pass_variants: int = 0
    filtered_variants: int = 0
    average_qual: float | None = None
    average_depth: float | None = None
    missing_genotypes: int = 0
    homozygous_count: int = 0
    heterozygous_count: int = 0
