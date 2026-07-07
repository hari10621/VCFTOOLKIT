"""Streaming VCF statistics."""

from __future__ import annotations

from .filters import classify_variant
from .models import VCFVariant, VariantStatistics
from .utils import first_numeric

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
BASES = {"A", "C", "G", "T"}


class StatisticsCollector:
    """Accumulates variant statistics without retaining variants."""

    def __init__(self) -> None:
        self.stats = VariantStatistics()
        self._qual_sum = 0.0
        self._qual_count = 0
        self._depth_sum = 0.0
        self._depth_count = 0

    def update(self, variant: VCFVariant) -> None:
        """Add a variant to the aggregate."""

        self.stats.total_variants += 1
        self.stats.chromosome_distribution[variant.chromosome] = (
            self.stats.chromosome_distribution.get(variant.chromosome, 0) + 1
        )
        if variant.filter == ["PASS"]:
            self.stats.pass_variants += 1
        else:
            self.stats.filtered_variants += 1

        variant_type = classify_variant(variant)
        if variant_type == "snp":
            self.stats.snp_count += 1
            self._update_snp_subtype(variant)
        elif variant_type == "insertion":
            self.stats.insertion_count += 1
        elif variant_type == "deletion":
            self.stats.deletion_count += 1
        elif variant_type == "structural":
            self.stats.structural_variant_count += 1

        if variant.qual is not None:
            self._qual_sum += variant.qual
            self._qual_count += 1
            self.stats.average_qual = self._qual_sum / self._qual_count

        depth = self._depth(variant)
        if depth is not None:
            self._depth_sum += depth
            self._depth_count += 1
            self.stats.average_depth = self._depth_sum / self._depth_count

        self._update_genotypes(variant)

    def snapshot(self) -> VariantStatistics:
        """Return the live statistics object."""

        return self.stats

    def _update_snp_subtype(self, variant: VCFVariant) -> None:
        ref = variant.reference.upper()
        for alt in variant.alternate:
            alt_upper = str(alt).upper()
            if ref not in BASES or alt_upper not in BASES:
                continue
            if (ref, alt_upper) in TRANSITIONS:
                self.stats.transition_count += 1
            else:
                self.stats.transversion_count += 1

    def _depth(self, variant: VCFVariant) -> float | None:
        depth = first_numeric(variant.info.get("DP"))
        if depth is not None:
            return depth
        sample_depths = [
            first_numeric(sample.get("DP"))
            for sample in variant.sample_values.values()
            if first_numeric(sample.get("DP")) is not None
        ]
        if not sample_depths:
            return None
        return sum(sample_depths) / len(sample_depths)

    def _update_genotypes(self, variant: VCFVariant) -> None:
        for sample in variant.sample_values.values():
            genotype = sample.get("GT")
            if genotype in {None, ".", "./.", ".|."}:
                self.stats.missing_genotypes += 1
                continue
            alleles = str(genotype).replace("|", "/").split("/")
            if "." in alleles:
                self.stats.missing_genotypes += 1
            elif len(set(alleles)) == 1:
                self.stats.homozygous_count += 1
            else:
                self.stats.heterozygous_count += 1
