import unittest
from vcf_reader.models import VCFVariant
from vcf_reader.inference import parse_ann_or_csq, extract_clinical_metadata, is_clinically_significant, generate_rule_based_report

class TestInference(unittest.TestCase):
    def test_parse_ann_snpeff(self) -> None:
        # Mock SnpEff variant
        variant = VCFVariant(
            chromosome="17",
            position=43044295,
            id="rs123456",
            reference="A",
            alternate=["G"],
            qual=100.0,
            filter=["PASS"],
            info={
                "ANN": "G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|NM_007294.3|Coding|11/23|c.98C>T|p.Ala33Val|182/1000|99/642|33/214||",
                "CLNSIG": "Pathogenic"
            }
        )
        
        parsed = parse_ann_or_csq(variant)
        self.assertEqual(parsed["gene"], "BRCA1")
        self.assertEqual(parsed["consequence"], "missense_variant")
        self.assertEqual(parsed["impact"], "MODERATE")
        self.assertEqual(parsed["hgvsc"], "c.98C>T")
        self.assertEqual(parsed["hgvsp"], "p.Ala33Val")
        
        metadata = extract_clinical_metadata(variant)
        self.assertEqual(metadata["gene"], "BRCA1")
        self.assertEqual(metadata["clnsig"], "Pathogenic")
        self.assertEqual(metadata["impact"], "MODERATE")
        
        self.assertTrue(is_clinically_significant(variant))

    def test_parse_csq_vep(self) -> None:
        # Mock VEP variant
        variant = VCFVariant(
            chromosome="7",
            position=117199646,
            id="rs789",
            reference="G",
            alternate=["A"],
            qual=99.0,
            filter=["PASS"],
            info={
                "CSQ": "A|stop_gained|HIGH|CFTR|ENSG00000001626|Transcript|ENST00000003084|Coding||||||||||||||||Pathogenic"
            }
        )
        
        parsed = parse_ann_or_csq(variant)
        self.assertEqual(parsed["gene"], "CFTR")
        self.assertEqual(parsed["consequence"], "stop_gained")
        self.assertEqual(parsed["impact"], "HIGH")
        
        metadata = extract_clinical_metadata(variant)
        self.assertEqual(metadata["gene"], "CFTR")
        self.assertEqual(metadata["clnsig"], "Pathogenic")
        self.assertEqual(metadata["impact"], "HIGH")
        
        self.assertTrue(is_clinically_significant(variant))

    def test_generate_report(self) -> None:
        variants_data = [
            {
                "chromosome": "17",
                "position": 43044295,
                "ref": "A",
                "alt": "G",
                "gene": "BRCA1",
                "id": "rs123456",
                "consequence": "missense_variant",
                "clnsig": "Pathogenic",
                "impact": "MODERATE",
                "hgvsc": "c.98C>T",
                "hgvsp": "p.Ala33Val",
                "disease": "Breast-ovarian cancer"
            }
        ]
        report = generate_rule_based_report(variants_data, "test.vcf")
        self.assertIn("GENOMIC INFERENCE REPORT", report)
        self.assertIn("BRCA1", report)
        self.assertIn("Breast-ovarian cancer", report)

if __name__ == "__main__":
    unittest.main()
