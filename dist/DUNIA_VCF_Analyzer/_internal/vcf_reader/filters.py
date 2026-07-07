"""Filtering and searching helpers."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any

from .models import VCFVariant, VariantFilter
from .utils import first_numeric


def classify_variant(variant: VCFVariant) -> str:
    """Classify a variant into a broad type."""

    if any(str(alt).startswith("<") and str(alt).endswith(">") for alt in variant.alternate):
        return "structural"
    ref_len = len(variant.reference)
    alt_lengths = [len(str(alt)) for alt in variant.alternate if alt]
    if ref_len == 1 and all(length == 1 for length in alt_lengths):
        return "snp"
    if any(length > ref_len for length in alt_lengths):
        return "insertion"
    if any(length < ref_len for length in alt_lengths):
        return "deletion"
    return "complex"


def variant_matches(variant: VCFVariant, criteria: VariantFilter) -> bool:
    """Return True when a variant satisfies all filter criteria."""

    if criteria.chromosome and variant.chromosome != criteria.chromosome:
        return False
    if criteria.position_start is not None and variant.position < criteria.position_start:
        return False
    if criteria.position_end is not None and variant.position > criteria.position_end:
        return False
    if criteria.min_qual is not None and (variant.qual is None or variant.qual < criteria.min_qual):
        return False
    if criteria.pass_only and variant.filter != ["PASS"]:
        return False
    if criteria.variant_type and classify_variant(variant) != criteria.variant_type.lower():
        return False
    if criteria.gene and not _contains_value(variant.info, criteria.gene):
        return False
    if criteria.min_depth is not None and _depth(variant) < criteria.min_depth:
        return False
    if criteria.min_allele_frequency is not None or criteria.max_allele_frequency is not None:
        af = _allele_frequency(variant)
        if af is None:
            return False
        if criteria.min_allele_frequency is not None and af < criteria.min_allele_frequency:
            return False
        if criteria.max_allele_frequency is not None and af > criteria.max_allele_frequency:
            return False
    for key, expected in criteria.info_fields.items():
        if key not in variant.info:
            return False
        if expected is not None and variant.info[key] != expected:
            return False
    if criteria.sample and criteria.sample not in variant.sample_values:
        return False
    return True


def filter_stream(variants: Iterable[VCFVariant], criteria: VariantFilter) -> Iterator[VCFVariant]:
    """Yield variants that match criteria."""

    for variant in variants:
        if variant_matches(variant, criteria):
            yield variant


def search_variant(variant: VCFVariant, **query: Any) -> bool:
    """Search a variant by common fields."""

    gene = query.get("gene") or query.get("gene_name")
    if gene and not _contains_value(variant.info, gene):
        return False
    chromosome = query.get("chromosome")
    if chromosome and variant.chromosome != chromosome:
        return False
    position = query.get("position")
    if position is not None and variant.position != int(position):
        return False
    variant_id = query.get("variant_id") or query.get("rsid") or query.get("rsID")
    if variant_id and variant.id != variant_id:
        return False
    reference = query.get("reference") or query.get("reference_allele")
    if reference and variant.reference != reference:
        return False
    alternate = query.get("alternate") or query.get("alternate_allele")
    if alternate and alternate not in variant.alternate:
        return False
    return True


def _contains_value(values: dict[str, Any], needle: str) -> bool:
    needle_lower = str(needle).lower()
    for value in values.values():
        if isinstance(value, (list, tuple)):
            if any(needle_lower == str(item).lower() for item in value):
                return True
        elif needle_lower == str(value).lower():
            return True
    return False


def _depth(variant: VCFVariant) -> int:
    info_depth = first_numeric(variant.info.get("DP"))
    if info_depth is not None:
        return int(info_depth)
    depths = [
        first_numeric(sample.get("DP"))
        for sample in variant.sample_values.values()
        if first_numeric(sample.get("DP")) is not None
    ]
    return int(sum(depths)) if depths else 0


def _allele_frequency(variant: VCFVariant) -> float | None:
    for key in ("AF", "FREQ"):
        value = first_numeric(variant.info.get(key))
        if value is not None:
            return value
    return None
