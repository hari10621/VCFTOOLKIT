"""Clinical inference and AI interpretation utilities for VCF variants."""

from __future__ import annotations

import json
import urllib.request
import urllib.error
from typing import Any
from .models import VCFVariant


def parse_ann_or_csq(variant: VCFVariant) -> dict[str, str]:
    """Parse SnpEff 'ANN' or Ensembl VEP 'CSQ' fields from variant INFO.
    
    Returns a dict with keys: gene, consequence, impact, hgvsc, hgvsp, clnsig
    """
    result = {
        "gene": "",
        "consequence": "",
        "impact": "",
        "hgvsc": "",
        "hgvsp": "",
        "clnsig": ""
    }
    
    # Try SnpEff ANN field first
    ann = variant.info.get("ANN")
    if ann:
        if not isinstance(ann, (list, tuple)):
            ann = [ann]
        for item in ann:
            parts = [p.strip() for p in str(item).split("|")]
            if len(parts) >= 11:
                # SnpEff format: Allele | Consequence | Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p
                if parts[3]: result["gene"] = parts[3]
                if parts[1]: result["consequence"] = parts[1]
                if parts[2]: result["impact"] = parts[2]
                if parts[9]: result["hgvsc"] = parts[9]
                if parts[10]: result["hgvsp"] = parts[10]
                break # Return first annotation
                
    # Try VEP CSQ field
    csq = variant.info.get("CSQ")
    if csq:
        if not isinstance(csq, (list, tuple)):
            csq = [csq]
        for item in csq:
            parts = [p.strip() for p in str(item).split("|")]
            # Standard VEP contains: Allele | Consequence | IMPACT | SYMBOL | Gene | Feature | ...
            if len(parts) >= 7:
                if parts[3]: result["gene"] = parts[3]
                if parts[1]: result["consequence"] = parts[1]
                if parts[2]: result["impact"] = parts[2]
                # Try to extract HGVSc and HGVSp if they are present in downstream positions
                # Often HGVSc is index 10 or 11, HGVSp is index 11 or 12 depending on VEP version
                for part in parts:
                    if part.startswith("c."):
                        result["hgvsc"] = part
                    elif part.startswith("p."):
                        result["hgvsp"] = part
                break
                
    return result


def extract_clinical_metadata(variant: VCFVariant) -> dict[str, Any]:
    """Extract clinical metadata from a VCF variant."""
    parsed = parse_ann_or_csq(variant)
    
    # Gene symbol
    gene = variant.info.get("GENE") or variant.info.get("Gene") or variant.info.get("SYMBOL") or parsed["gene"] or ""
    if isinstance(gene, (list, tuple)) and gene:
        gene = str(gene[0])
    else:
        gene = str(gene)
        
    # Clinical significance (ClinVar)
    clnsig = variant.info.get("CLNSIG") or variant.info.get("CLN_SIG") or variant.info.get("clinical_significance") or parsed["clnsig"]
    if clnsig:
        if isinstance(clnsig, (list, tuple)):
            clnsig_str = "/".join(str(x) for x in clnsig)
        else:
            clnsig_str = str(clnsig)
    else:
        clnsig_str = ""
        
    # Try to extract clinical significance from other fields or ANN/CSQ details
    if not clnsig_str:
        # Check if any field in info has pathogenic, benign, VUS etc
        for k, v in variant.info.items():
            val_str = str(v).lower()
            if any(term in val_str for term in ["pathogenic", "likely_pathogenic", "uncertain_significance", "vus", "benign"]):
                if k.startswith("CLN"):
                    clnsig_str = str(v)
                    break

    # If still empty, check ANN/CSQ strings for clinical significance terms
    if not clnsig_str:
        ann_csq = variant.info.get("ANN") or variant.info.get("CSQ")
        if ann_csq:
            if not isinstance(ann_csq, (list, tuple)):
                ann_csq = [ann_csq]
            for item in ann_csq:
                for part in str(item).split("|"):
                    part_lower = part.lower()
                    if "pathogenic" in part_lower:
                        clnsig_str = "Pathogenic"
                        break
                    elif "likely_pathogenic" in part_lower:
                        clnsig_str = "Likely pathogenic"
                        break
                    elif "uncertain_significance" in part_lower or "vus" in part_lower:
                        clnsig_str = "Uncertain significance"
                        break

    # Consequence and Impact
    impact = variant.info.get("IMPACT") or variant.info.get("impact") or parsed["impact"] or ""
    consequence = variant.info.get("CONSEQUENCE") or variant.info.get("Consequence") or variant.info.get("MC") or parsed["consequence"] or ""
    
    impact_str = str(impact) if impact else ""
    conseq_str = str(consequence) if consequence else ""
    
    if not impact_str:
        # Infer impact from consequence terms
        c_lower = conseq_str.lower()
        if any(term in c_lower for term in ["stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost"]):
            impact_str = "HIGH"
        elif any(term in c_lower for term in ["missense", "inframe", "protein_altering"]):
            impact_str = "MODERATE"
        elif any(term in c_lower for term in ["synonymous", "intron", "splice_region"]):
            impact_str = "LOW"
            
    # Disease Name / Phenotype
    disease = variant.info.get("CLNDN") or variant.info.get("CLNDISDB") or ""
    if isinstance(disease, (list, tuple)) and disease:
        disease = "/".join(str(d) for d in disease)
    else:
        disease = str(disease)
        
    return {
        "gene": gene or "Unknown",
        "clnsig": clnsig_str or "Not annotated",
        "impact": impact_str or "MODIFIER",
        "consequence": conseq_str or "Sequence variant",
        "disease": disease or "Not specified",
        "hgvsc": parsed["hgvsc"] or f"c.{variant.position}{variant.reference}>{','.join(variant.alternate)}",
        "hgvsp": parsed["hgvsp"] or ""
    }


def is_clinically_significant(variant: VCFVariant) -> bool:
    """Determine if a variant is potentially clinically significant."""
    metadata = extract_clinical_metadata(variant)
    clnsig_lower = metadata["clnsig"].lower()
    impact_upper = metadata["impact"].upper()
    
    # Pathogenic or likely pathogenic in ClinVar
    if "pathogenic" in clnsig_lower:
        return True
    # High impact functional consequence (stop gained, frameshift, splice site disruption)
    if impact_upper == "HIGH":
        return True
    # Non-synonymous missense variants in key genes or with high quality
    if impact_upper == "MODERATE" and variant.qual and variant.qual >= 50:
        return True
        
    return False


def call_gemini_api(api_key: str, model: str, prompt: str) -> str:
    """Query the Google Gemini API using native urllib.
    
    This avoids dependencies on external SDKs in the compiled executable.
    """
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={api_key}"
    headers = {"Content-Type": "application/json"}
    
    payload = {
        "contents": [
            {
                "parts": [
                    {"text": prompt}
                ]
            }
        ],
        "generationConfig": {
            "temperature": 0.2,
            "maxOutputTokens": 2048
        }
    }
    
    try:
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(url, data=data, headers=headers, method="POST")
        
        with urllib.request.urlopen(req, timeout=30) as response:
            res_data = json.loads(response.read().decode("utf-8"))
            
            # Extract text from the candidate response
            candidates = res_data.get("candidates", [])
            if candidates:
                parts = candidates[0].get("content", {}).get("parts", [])
                if parts:
                    return str(parts[0].get("text", ""))
            
            return f"Error: Received empty response or unexpected format from Gemini. Full response: {res_data}"
            
    except urllib.error.HTTPError as exc:
        try:
            error_body = exc.read().decode("utf-8")
            error_json = json.loads(error_body)
            error_msg = error_json.get("error", {}).get("message", str(exc))
            return f"Gemini API HTTP Error: {exc.code} - {error_msg}"
        except Exception:
            return f"Gemini API HTTP Error: {exc}"
    except Exception as exc:
        return f"Gemini API Request Failed: {exc}"


def generate_rule_based_report(variants_data: list[dict[str, Any]], filename: str) -> str:
    """Generate a clean markdown clinical report offline using heuristic rules."""
    
    pathogenic_count = sum(1 for v in variants_data if "pathogenic" in v["clnsig"].lower() and "likely" not in v["clnsig"].lower())
    likely_pathogenic_count = sum(1 for v in variants_data if "likely pathogenic" in v["clnsig"].lower())
    high_impact_count = sum(1 for v in variants_data if v["impact"].upper() == "HIGH")
    vus_count = sum(1 for v in variants_data if "uncertain" in v["clnsig"].lower() or "vus" in v["clnsig"].lower())
    
    md = []
    md.append(f"# GENOMIC INFERENCE REPORT (OFFLINE)")
    md.append(f"**Source File:** `{filename}`")
    md.append(f"**Analysis Mode:** Rule-Based Annotation Parsing (Offline)")
    md.append("")
    md.append("## Executive Summary")
    md.append(f"A genomic scan was performed on the variant call dataset. A total of **{len(variants_data)}** high-priority variants were extracted based on ClinVar pathogenicity designations and high functional consequence impact scores.")
    md.append("")
    md.append(f"- **Pathogenic Variants:** {pathogenic_count}")
    md.append(f"- **Likely Pathogenic Variants:** {likely_pathogenic_count}")
    md.append(f"- **High-Impact Consequence Variants (Frameshifts, Nonsense):** {high_impact_count}")
    md.append(f"- **Variants of Uncertain Significance (VUS):** {vus_count}")
    md.append("")
    
    if pathogenic_count + likely_pathogenic_count > 0:
        md.append("> [!WARNING]")
        md.append(f"> **Clinical Findings Detected:** This sample contains {pathogenic_count + likely_pathogenic_count} variant(s) annotated as Clinically Significant (Pathogenic/Likely Pathogenic) in ClinVar. Review details below.")
        md.append("")
    else:
        md.append("> [!NOTE]")
        md.append("> No clear pathogenic or likely pathogenic variants were detected under current filters. However, other variants of interest are listed below.")
        md.append("")

    md.append("## High-Priority Variants Table")
    md.append("| Chromosome | Position | Gene | ID | Consequence | ClinVar Significance | Impact |")
    md.append("|---|---|---|---|---|---|---|")
    for v in variants_data[:50]: # Limit to 50 in table for readability
        md.append(f"| {v['chromosome']} | {v['position']} | **{v['gene']}** | {v['id']} | {v['consequence']} | {v['clnsig']} | `{v['impact']}` |")
    md.append("")
    
    md.append("## Gene-Level Details and Interpretations")
    seen_genes = set()
    for v in variants_data:
        gene = v["gene"]
        if gene == "Unknown" or gene in seen_genes:
            continue
        seen_genes.add(gene)
        md.append(f"### Gene: {gene}")
        md.append(f"Variants in this gene can be associated with phenotypes including: *{v['disease']}*.")
        md.append(f"- **Representative Variant:** `{v['chromosome']}:{v['position']} {v['ref']}>{v['alt']}`")
        md.append(f"- **HGVS nomenclature:** `{v['hgvsc']}` {f'({v['hgvsp']})' if v['hgvsp'] else ''}")
        md.append(f"- **Functional Impact:** {v['consequence']} (Impact level: `{v['impact']}`)")
        md.append(f"- **ClinVar Status:** **{v['clnsig']}**")
        md.append("")
        
    md.append("---")
    md.append("## Disclaimer")
    md.append("This report is generated automatically by the DUNIA VCF Analyzer based on annotations found within the VCF file. It is for research purposes only and does not constitute formal medical or clinical advice. All clinical findings should be validated by an accredited diagnostic laboratory and interpreted by a certified clinical geneticist.")
    
    return "\n".join(md)


def build_gemini_prompt(variants_data: list[dict[str, Any]], filename: str, report_type: str) -> str:
    """Build a detailed prompt for the Gemini AI clinical interpreter."""
    
    # Format the top variants into JSON for the prompt
    formatted_variants = []
    for v in variants_data[:30]: # Send top 30 significant variants to fit token limits comfortably
        formatted_variants.append({
            "Chr": v["chromosome"],
            "Pos": v["position"],
            "Ref": v["ref"],
            "Alt": v["alt"],
            "Gene": v["gene"],
            "rsID": v["id"],
            "Consequence": v["consequence"],
            "ClinVarSig": v["clnsig"],
            "Impact": v["impact"],
            "HGVS_c": v["hgvsc"],
            "HGVS_p": v["hgvsp"],
            "DiseaseName": v["disease"]
        })
        
    variants_json_str = json.dumps(formatted_variants, indent=2)
    
    prompt = f"""You are a professional clinical geneticist and expert medical bioinformatics AI.
You have been provided with a list of clinically significant or high-impact genetic variants parsed from a patient's VCF file (filename: {filename}).

Here is the parsed variants data in JSON format:
{variants_json_str}

Please generate a professional, structured clinical interpretation report based on these variants.
Report Type Requested: {report_type}

Follow these structural guidelines:
1. Start with a professional header: "GENOMIC CLINICAL INTERPRETATION REPORT".
2. **Executive Summary**: Synthesize the key findings. Highlight if there are any Pathogenic or Likely Pathogenic mutations (especially in well-known genes like BRCA1/2, TP53, CFTR, etc.) and what diseases/syndromes they relate to.
3. **Primary Findings Table**: A neat markdown table of the most critical variants.
4. **Variant-by-Variant Clinical Analysis**: Write a detailed clinical interpretation for the most critical variants. For each critical variant:
   - State the Gene, HGVS nomenclature (c. and p. if available), and rsID.
   - Describe the biological mechanism of the variant (e.g. missense, stop gained leading to nonsense-mediated decay).
   - Explain the clinical significance, citing ClinVar classifications.
   - Mention the associated disease/phenotype and pattern of inheritance (e.g., autosomal dominant, autosomal recessive) if known.
5. **Next Steps and Recommendations**: Suggest appropriate clinical follow-ups (e.g., genetic counseling, confirmatory Sanger sequencing, family member testing, specific clinical screenings).
6. **Disclaimer**: Include a standard medical disclaimer stating that this report is generated by an AI assistant for research/educational purposes and must be verified by a licensed healthcare professional.

Use clean, professional Markdown formatting. Do not output HTML tags. Keep the tone clinical, precise, and objective.
"""
    return prompt
