from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from vcf_reader import VCFReader, VCFValidationError, VariantFilter
from vcf_reader.models import ParserBackend


VCF_TEXT = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\trs1\tA\tG\t50\tPASS\tDP=10;AF=0.5;GENE=BRCA1\tGT:DP\t0/1:5\t1/1:5
1\t120\t.\tAT\tA\t10\tLowQual\tDP=6;AF=0.1;GENE=TP53\tGT:DP\t0/0:6\t./.:.
2\t200\trs2\tC\tT\t90\tPASS\tDP=12;AF=0.2;GENE=BRCA2\tGT:DP\t0/1:7\t0/0:5
"""


class VCFReaderTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.tmp.name)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def write_vcf(self, text: str = VCF_TEXT, name: str = "sample.vcf") -> Path:
        path = self.tmp_path / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_read_header_and_stream_small_vcf(self) -> None:
        reader = VCFReader(self.write_vcf(), backend=ParserBackend.MANUAL)

        self.assertTrue(reader.validate())
        header = reader.read_header()
        self.assertEqual(header.version, "VCFv4.2")
        self.assertEqual(header.reference, "GRCh38")
        self.assertIn("DP", header.info_definitions)
        self.assertEqual(header.samples, ["S1", "S2"])

        variants = list(reader.stream_variants())
        self.assertEqual(len(variants), 3)
        self.assertEqual(variants[0].chromosome, "1")
        self.assertEqual(variants[0].position, 100)
        self.assertEqual(variants[0].alternate, ["G"])
        self.assertEqual(variants[0].sample_values["S1"]["GT"], "0/1")

    def test_statistics_are_collected_while_streaming(self) -> None:
        reader = VCFReader(self.write_vcf(), backend=ParserBackend.MANUAL)

        stats = reader.get_statistics(refresh=True)

        self.assertEqual(stats.total_variants, 3)
        self.assertEqual(stats.snp_count, 2)
        self.assertEqual(stats.deletion_count, 1)
        self.assertEqual(stats.transition_count, 2)
        self.assertEqual(stats.pass_variants, 2)
        self.assertEqual(stats.filtered_variants, 1)
        self.assertEqual(stats.chromosome_distribution, {"1": 2, "2": 1})
        self.assertEqual(stats.missing_genotypes, 1)
        self.assertEqual(stats.homozygous_count, 3)
        self.assertEqual(stats.heterozygous_count, 2)
        self.assertAlmostEqual(stats.average_depth or 0, 28 / 3)

    def test_filtering_and_searching(self) -> None:
        reader = VCFReader(self.write_vcf(), backend=ParserBackend.MANUAL)

        filtered = list(reader.filter_variants(VariantFilter(chromosome="1", min_qual=20, pass_only=True)))
        self.assertEqual([variant.id for variant in filtered], ["rs1"])

        search_results = list(reader.search(gene="BRCA2", limit=1))
        self.assertEqual(len(search_results), 1)
        self.assertEqual(search_results[0].id, "rs2")

    def test_export_csv_json_tsv_and_vcf(self) -> None:
        reader = VCFReader(self.write_vcf(), backend=ParserBackend.MANUAL)
        reader.read_header()

        csv_path = self.tmp_path / "out.csv"
        json_path = self.tmp_path / "out.json"
        tsv_path = self.tmp_path / "out.tsv"
        vcf_path = self.tmp_path / "out.vcf"

        self.assertEqual(reader.export(csv_path, "csv"), 3)
        self.assertEqual(reader.export(json_path, "json"), 3)
        self.assertEqual(reader.export(tsv_path, "tsv"), 3)
        self.assertEqual(reader.export(vcf_path, "vcf", VariantFilter(pass_only=True)), 2)

        self.assertIn("chromosome,position", csv_path.read_text(encoding="utf-8"))
        self.assertIn('"chromosome": "1"', json_path.read_text(encoding="utf-8"))
        self.assertIn("#CHROM\tPOS\tID", vcf_path.read_text(encoding="utf-8"))

    def test_compressed_vcf(self) -> None:
        path = self.tmp_path / "sample.vcf.gz"
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(VCF_TEXT)

        reader = VCFReader(path, backend=ParserBackend.MANUAL)

        self.assertTrue(reader.validate())
        self.assertEqual(len(list(reader.stream_variants())), 3)

    def test_corrupted_vcf_raises_validation_error(self) -> None:
        path = self.write_vcf("##fileformat=VCFv4.2\n#CHROM\tPOS\n1\tbad\n", "bad.vcf")
        reader = VCFReader(path, backend=ParserBackend.MANUAL)

        with self.assertRaises(VCFValidationError):
            reader.validate()

    def test_missing_header_raises_validation_error(self) -> None:
        path = self.write_vcf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", "missing.vcf")
        reader = VCFReader(path, backend=ParserBackend.MANUAL)

        with self.assertRaisesRegex(VCFValidationError, "Missing required ##fileformat"):
            reader.validate()

    def test_large_vcf_is_streamed(self) -> None:
        records = [
            f"1\t{i}\trs{i}\tA\tG\t50\tPASS\tDP=8;AF=0.1;GENE=GENE{i % 3}\tGT:DP\t0/1:8\t0/0:8"
            for i in range(1, 5001)
        ]
        text = VCF_TEXT.split("#CHROM", 1)[0]
        text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        text += "\n".join(records) + "\n"
        reader = VCFReader(self.write_vcf(text), backend=ParserBackend.MANUAL)

        count = 0
        for _variant in reader.stream_variants():
            count += 1

        self.assertEqual(count, 5000)
        self.assertEqual(reader.progress().current_variant_count, 5000)


if __name__ == "__main__":
    unittest.main()
