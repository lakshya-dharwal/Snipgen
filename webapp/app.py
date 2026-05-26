"""SnipGen FastAPI web application — async job queue edition."""

import json
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, UploadFile, HTTPException, Query, Body
from fastapi.responses import HTMLResponse, JSONResponse

from snipgen.analysis.base_editor import analyze_base_editing
from snipgen.analysis.cloning_primers import design_cloning_oligos, design_all_vectors
from snipgen.analysis.isoform_analyzer import analyze_guide_isoforms
from snipgen.filters.gnomad_filter import check_guide_gnomad
from snipgen.filters.pam_filter import PAM_REGISTRY
from snipgen.pipeline import PipelineConfig, SnipGenPipeline
from snipgen.scoring.clinvar_annotator import annotate_gene as clinvar_annotate_gene
from webapp.crispor_client import (
    submit_sequence as crispor_submit,
    fetch_scores as crispor_fetch,
    crispor_to_offtarget_score,
)
from webapp.job_queue import queue, JobStatus

app = FastAPI(title="SnipGen", description="AI-driven CRISPR guide RNA design")

_static = Path(__file__).resolve().parent / "static"
_ENTREZ_EMAIL = "snipgen-tool@noreply.asu.edu"

_GENE_ACCESSIONS: dict[str, dict[str, str]] = {
    "human": {
        "TP53":  "NM_000546", "BRCA1": "NM_007294", "BRCA2": "NM_000059",
        "EGFR":  "NM_005228", "KRAS":  "NM_004985", "PTEN":  "NM_000314",
        "MYC":   "NM_002467", "VEGFA": "NM_001171627", "PCSK9": "NM_174936",
        "HBB":   "NM_000518", "DMD":   "NM_004006",  "CFTR":  "NM_000492",
        "APOE":  "NM_000041", "ACE2":  "NM_021804",  "STAT3": "NM_139276",
    },
    "mouse": {
        "Trp53": "NM_011640", "Brca1": "NM_009764", "Kras": "NM_021284",
        "Egfr":  "NM_207655", "Myc":   "NM_010849", "Hbb":  "NM_008220",
        "Dmd":   "NM_007868",
    },
}

ORGANISM_TAXIDS = {
    "human": "9606", "mouse": "10090", "zebrafish": "7955",
    "rat": "10116", "drosophila": "7227", "c_elegans": "6239",
}


# ── Helpers ────────────────────────────────────────────────────────────────────

def _fetch_sequence_entrez(gene: str, organism: str) -> tuple[str, str]:
    from Bio import Entrez
    Entrez.email = _ENTREZ_EMAIL
    Entrez.tool  = "snipgen"
    org_lower    = organism.lower()
    gene_upper   = gene.upper()
    accession: Optional[str] = None

    if org_lower in _GENE_ACCESSIONS and gene_upper in _GENE_ACCESSIONS[org_lower]:
        accession = _GENE_ACCESSIONS[org_lower][gene_upper]

    if accession is None:
        taxid = ORGANISM_TAXIDS.get(org_lower, "9606")
        term = f"{gene}[Gene Name] AND {taxid}[Taxonomy ID] AND mRNA[Filter] AND RefSeq[Filter]"
        try:
            h = Entrez.esearch(db="nuccore", term=term, retmax=5)
            rec = Entrez.read(h); h.close()
        except Exception as exc:
            raise ValueError(f"NCBI search failed: {exc}")
        ids = rec.get("IdList", [])
        if not ids:
            raise ValueError(f"No RefSeq mRNA found for '{gene}' in {organism}.")
        accession = ids[0]

    try:
        h = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
        fasta = h.read(); h.close()
    except Exception as exc:
        raise ValueError(f"NCBI fetch failed for {accession}: {exc}")

    if not fasta.strip():
        raise ValueError(f"Empty sequence for {accession}.")
    return fasta, accession


def _compute_decision_margin(candidates: list) -> dict:
    """
    Compute confidence margin between top guide and alternatives.
    Inspired by ADMET confidence/margin thinking: a score alone is insufficient
    without knowing how far it is from the next-best candidate.
    """
    if not candidates:
        return {}
    if len(candidates) == 1:
        return {
            "margin": 100.0,
            "label": "SINGLE",
            "top_score": candidates[0].final_score,
            "second_score": None,
            "third_score": None,
            "advice": "Only one guide passed all filters. Experimental validation is essential.",
        }

    top    = candidates[0].final_score
    second = candidates[1].final_score
    third  = candidates[2].final_score if len(candidates) > 2 else None
    margin = round(top - second, 1)

    if margin >= 15:
        label  = "STRONG"
        advice = (
            f"Top guide scores {margin} pts above the next candidate "
            f"({top} vs {second}). Ranking is reliable — high confidence in top selection."
        )
    elif margin >= 7:
        label  = "MODERATE"
        advice = (
            f"Top two guides are {margin} pts apart ({top} vs {second}). "
            "Ranking is moderately confident. Consider validating the top 2 guides wet-lab."
        )
    else:
        label  = "WEAK"
        advice = (
            f"Top guides are within {margin} pts of each other ({top} vs {second}). "
            "Scores are statistically close — recommend wet-lab validation of top 3–5 guides "
            "before committing to a single candidate."
        )

    return {
        "margin": margin,
        "label": label,
        "top_score": top,
        "second_score": second,
        "third_score": third,
        "advice": advice,
    }


def _build_explainability(candidate, rank: int) -> list[str]:
    """
    Return human-readable bullet points explaining why this guide was ranked here.
    Each bullet answers: what factor, what value, what it means.
    """
    bd     = candidate.score_breakdown or {}
    seq    = candidate.sequence.upper()
    gc_pct = round(candidate.gc_content * 100, 1)
    points: list[str] = []

    # On-target score
    ots = candidate.on_target_score
    if ots >= 70:
        points.append(
            f"✅ On-target efficiency {ots}/100 — high confidence from Azimuth-class model "
            f"(GBR trained on 5,310 real CRISPR screen guides, Spearman r=0.556)"
        )
    elif ots >= 45:
        points.append(
            f"⚠ On-target efficiency {ots}/100 — moderate; consider experimental validation "
            f"(T7E1 or TIDE assay) before committing"
        )
    else:
        points.append(
            f"✗ On-target efficiency {ots}/100 — lower than alternatives; ranked here because "
            f"other factors compensate, but treat with caution"
        )

    # GC content
    if 45 <= gc_pct <= 65:
        points.append(f"✅ GC content {gc_pct}% — in optimal range (45–65%) for stable Cas9 binding")
    elif 40 <= gc_pct <= 70:
        points.append(f"⚠ GC content {gc_pct}% — acceptable but slightly outside peak 45–65%")
    else:
        points.append(
            f"✗ GC content {gc_pct}% — outside optimal 40–70% range, may reduce efficiency"
        )

    # Isoform coverage
    iso = bd.get("isoform") or getattr(candidate, "isoform", None)
    if iso and iso.get("isoform_checked"):
        lbl   = iso.get("label", "")
        n_hit = iso.get("n_hits", 0)
        n_tot = iso.get("n_total", 0)
        if lbl in ("PAN-ISOFORM", "BROAD"):
            points.append(
                f"✅ Isoform coverage: {lbl} — found in {n_hit}/{n_tot} RefSeq transcripts; "
                f"ensures complete gene-level knockout"
            )
        elif lbl == "SELECTIVE":
            points.append(
                f"⚠ Isoform coverage: SELECTIVE ({n_hit}/{n_tot} transcripts) — "
                f"verify this covers your target isoform"
            )
        else:
            points.append(
                f"✗ Isoform coverage: SINGLE-ISOFORM (1/{n_tot} transcripts) — "
                f"incomplete knockdown; may miss disease-relevant splice forms"
            )

    # gnomAD seed-SNP check
    gnomad = bd.get("gnomad") or getattr(candidate, "gnomad", None)
    if gnomad and gnomad.get("gnomad_checked"):
        risk = gnomad.get("risk_level", "NONE")
        af   = gnomad.get("af_max", 0)
        if risk == "NONE":
            points.append(
                "✅ No common population variants in seed region (gnomAD v4) — "
                "sequence stable across human populations"
            )
        elif risk in ("LOW", "MODERATE"):
            points.append(
                f"⚠ Minor population variant overlap (gnomAD AF ≈ {af:.3f}, risk: {risk}) — "
                f"monitor in diverse cell lines"
            )
        else:
            points.append(
                f"✗ Common seed-region SNP detected (gnomAD AF {af:.3f}) — "
                f"Cas9 binding may be disrupted in a fraction of samples"
            )

    # Off-target risk
    ot_score = candidate.off_target_score
    if ot_score >= 75:
        points.append(f"✅ Low predicted off-target risk (heuristic score {ot_score}/100)")
    elif ot_score >= 50:
        points.append(
            f"⚠ Moderate off-target risk (score {ot_score}/100) — "
            f"await genome-wide CRISPOR results before use"
        )
    else:
        points.append(
            f"✗ Higher off-target risk (score {ot_score}/100) — "
            f"validate with GUIDE-seq or NGS before experimental use"
        )

    # Rank label
    if rank == 1:
        points.append("🏆 Highest composite score across all evaluated criteria in this run")
    elif rank <= 3:
        points.append(f"📊 Ranked #{rank} — strong alternative; include in experimental shortlist")

    return points


def _compute_do_not_use(candidate) -> dict:
    """
    Apply strict red-flag logic. A guide that fails ANY of these criteria
    should not be used without careful review, regardless of composite score.
    """
    bd      = candidate.score_breakdown or {}
    seq     = candidate.sequence.upper()
    gc      = candidate.gc_content
    reasons: list[str] = []

    # GC extremes
    if gc > 0.75:
        reasons.append(
            f"GC content {gc*100:.0f}% is extremely high (>75%) — "
            f"guides form secondary structures, reducing Cas9 loading efficiency"
        )
    if gc < 0.35:
        reasons.append(
            f"GC content {gc*100:.0f}% is extremely low (<35%) — "
            f"insufficient melting temperature for stable R-loop formation"
        )

    # Homopolymer run ≥5
    for base in "ACGT":
        if base * 5 in seq:
            reasons.append(
                f"Homopolymer run detected ({base}×5+) — "
                f"severely reduces Pol III transcription efficiency and Cas9 cleavage"
            )
            break

    # Seed-region SNP
    gnomad = bd.get("gnomad") or getattr(candidate, "gnomad", None)
    if gnomad and gnomad.get("risk_level") == "HIGH":
        af = gnomad.get("af_max", 0)
        reasons.append(
            f"Common population SNP in seed region (gnomAD AF {af:.3f}) — "
            f"Cas9 binding will be disrupted in a significant fraction of samples"
        )

    # Single-isoform with many alternatives
    iso = bd.get("isoform") or getattr(candidate, "isoform", None)
    if iso and iso.get("label") == "SINGLE-ISOFORM" and (iso.get("n_total", 0) or 0) >= 4:
        reasons.append(
            f"Targets only 1 of {iso['n_total']} isoforms — "
            f"incomplete knockdown; gene products from other isoforms unaffected"
        )

    # Safety tier AVOID
    if candidate.safety_label == "AVOID":
        reasons.append(
            "Safety tier AVOID — predicted off-target profile exceeds acceptable threshold"
        )

    # PolyT (U6 terminator signal)
    if "TTTTT" in seq:
        reasons.append(
            "Five consecutive Ts detected — acts as RNA Pol III termination signal, "
            "truncating guide RNA transcription"
        )

    return {"do_not_use": len(reasons) > 0, "reasons": reasons}


def _build_rejected_examples(rejected_candidates: list, n: int = 3) -> list[dict]:
    """
    Expose the top-N rejected guides with clear rejection reasons.
    These are intentionally surfaced to show the system's filter logic transparently.
    """
    examples = []
    for c in rejected_candidates[:n]:
        bd      = c.score_breakdown or {}
        seq     = c.sequence.upper()
        gc_pct  = round(c.gc_content * 100, 1)
        reasons: list[str] = []

        if c.gc_content > 0.70:
            reasons.append(f"GC content {gc_pct}% exceeds 70% maximum")
        elif c.gc_content < 0.40:
            reasons.append(f"GC content {gc_pct}% below 40% minimum")

        for base in "ACGT":
            if base * 5 in seq:
                reasons.append(f"Homopolymer run ({base}×5) detected")
                break

        if "TTTTT" in seq:
            reasons.append("Poly-T run → Pol III termination signal")

        seed_gc = bd.get("seed_gc", None)
        if seed_gc is not None and seed_gc > 0.75:
            reasons.append(f"Seed-region GC {seed_gc*100:.0f}% too high (>75%)")

        flt = bd.get("filter_reason", "")
        if flt and flt not in " ".join(reasons):
            reasons.append(flt)

        if not reasons:
            reasons.append("Did not meet filter thresholds (GC, PAM, seed quality)")

        examples.append({
            "sequence":          c.sequence,
            "pam":               c.pam,
            "gc_content":        gc_pct,
            "safety_label":      c.safety_label,
            "final_score":       c.final_score,
            "rejection_reasons": reasons,
            "decision":          "REJECTED",
        })
    return examples


def _run_pipeline(fasta_bytes: bytes, cas_variant: str, guide_length: int,
                  min_gc: float, max_gc: float, top_n: int,
                  organism: str, gene_symbol: str = "") -> dict:
    """
    Run the full SnipGen pipeline. Executed in a background thread by job_queue.
    Returns the JSON-serialisable result dict.
    """
    fasta_text = fasta_bytes.decode("utf-8", errors="replace")

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="wb") as tmp:
        tmp.write(fasta_bytes)
        tmp_path = Path(tmp.name)

    crispor_batch_id: Optional[str] = None

    try:
        with tempfile.TemporaryDirectory() as out_dir:
            config = PipelineConfig(
                fasta_path=tmp_path,
                output_dir=out_dir,
                output_formats=["json"],
                cas_variant=cas_variant,
                guide_length=guide_length,
                min_gc=min_gc,
                max_gc=max_gc,
                top_n=top_n,
            )
            pipeline = SnipGenPipeline(config)
            pipeline.run()
            output = json.loads((Path(out_dir) / "candidates.json").read_text())

        gene = (gene_symbol or "").strip().upper()
        candidates = output.get("candidates", [])

        # ── Base editing analysis (always run — pure sequence analysis) ───────
        for c in candidates:
            seq    = c.get("sequence", "")
            strand = c.get("strand", "+")
            try:
                c["base_edit"] = analyze_base_editing(seq, strand)
            except Exception:
                pass

        # ── Cloning primers (always run — pure sequence, no network) ──────────
        for c in candidates:
            seq = c.get("sequence", "")
            try:
                c["cloning"] = design_all_vectors(seq)
            except Exception:
                pass

        # ── ClinVar gene annotation (gene-search mode) ────────────────────────
        if gene:
            try:
                c_ann = clinvar_annotate_gene(gene)
                output["clinvar_gene"] = c_ann
            except Exception:
                output["clinvar_gene"] = {}

        # ── gnomAD + isoform annotation (gene-search mode, human only) ────────
        if gene and organism.lower() == "human":
            for c in candidates:
                seq    = c.get("sequence", "")
                start  = c.get("start", 0)
                end    = c.get("end", start + len(seq))
                strand = c.get("strand", "+")

                try:
                    c["gnomad"] = check_guide_gnomad(
                        guide_seq=seq,
                        gene_symbol=gene,
                        guide_start_in_mrna=start,
                        guide_end_in_mrna=end,
                        strand=strand,
                        organism=organism,
                    )
                except Exception:
                    c["gnomad"] = {"gnomad_checked": False, "flag": "gnomAD check failed"}

            # Isoform — cap at top 8 guides (transcript FASTA cached after first)
            for c in candidates[:8]:
                seq = c.get("sequence", "")
                try:
                    c["isoform"] = analyze_guide_isoforms(
                        guide_seq=seq,
                        gene_symbol=gene,
                        organism=organism,
                    )
                except Exception:
                    c["isoform"] = {"isoform_checked": False, "flag": "Isoform check failed"}

        # ── CRISPOR submission (fire-and-forget) ──────────────────────────────
        try:
            crispor_batch_id = crispor_submit(fasta_text, organism, cas_variant)
        except Exception:
            crispor_batch_id = None

        output["crispor_batch_id"] = crispor_batch_id
        output["crispor_genome"]   = organism
        output["gene_symbol"]      = gene

        # ── Decision margin ────────────────────────────────────────────────
        # We need GRNACandidate objects for margin computation; re-run pipeline
        # to get them — but candidates list is already dicts from JSON.
        # Compute margin from the dict list directly.
        cand_dicts  = output.get("candidates", [])
        scores      = [c.get("final_score", 0) for c in cand_dicts]
        if len(scores) >= 2:
            margin_val = round(scores[0] - scores[1], 1)
            if margin_val >= 15:
                margin_label  = "STRONG"
                margin_advice = (
                    f"Top guide scores {margin_val} pts above the next candidate "
                    f"({scores[0]} vs {scores[1]}). Ranking is reliable — high confidence in top selection."
                )
            elif margin_val >= 7:
                margin_label  = "MODERATE"
                margin_advice = (
                    f"Top two guides are {margin_val} pts apart ({scores[0]} vs {scores[1]}). "
                    "Consider validating the top 2 guides experimentally."
                )
            else:
                margin_label  = "WEAK"
                margin_advice = (
                    f"Top guides are within {margin_val} pts of each other ({scores[0]} vs {scores[1]}). "
                    "Scores are statistically close — recommend wet-lab validation of top 3–5 guides."
                )
        elif len(scores) == 1:
            margin_val    = 100.0
            margin_label  = "SINGLE"
            margin_advice = "Only one guide passed all filters. Experimental validation is essential."
        else:
            margin_val, margin_label, margin_advice = 0.0, "NONE", ""

        output["decision_margin"] = {
            "margin":       margin_val,
            "label":        margin_label,
            "advice":       margin_advice,
            "top_score":    scores[0] if scores else None,
            "second_score": scores[1] if len(scores) > 1 else None,
            "third_score":  scores[2] if len(scores) > 2 else None,
        }

        # ── Per-candidate explainability + do-not-use flags ───────────────
        for rank_i, c in enumerate(cand_dicts):
            bd      = c.get("score_breakdown") or {}
            seq     = c.get("sequence", "").upper()
            gc      = c.get("gc_content", 0.5)
            gc_pct  = round(gc * 100, 1)
            ots     = c.get("on_target_score", 50)
            ot_sc   = c.get("off_target_score", 50)

            # Explainability bullets
            pts: list[str] = []
            if ots >= 70:
                pts.append(f"✅ On-target efficiency {ots}/100 — high confidence (Azimuth ML, Spearman 0.556)")
            elif ots >= 45:
                pts.append(f"⚠ On-target efficiency {ots}/100 — moderate; validate experimentally")
            else:
                pts.append(f"✗ On-target efficiency {ots}/100 — lower than alternatives")

            if 45 <= gc_pct <= 65:
                pts.append(f"✅ GC content {gc_pct}% — optimal range (45–65%) for stable Cas9 binding")
            elif 40 <= gc_pct <= 70:
                pts.append(f"⚠ GC content {gc_pct}% — acceptable but slightly outside peak 45–65%")
            else:
                pts.append(f"✗ GC content {gc_pct}% — outside optimal 40–70% range")

            iso = bd.get("isoform") or c.get("isoform")
            if iso and iso.get("isoform_checked"):
                lbl  = iso.get("label", "")
                n_h  = iso.get("n_hits", 0)
                n_t  = iso.get("n_total", 0)
                if lbl in ("PAN-ISOFORM", "BROAD"):
                    pts.append(f"✅ Isoform coverage: {lbl} ({n_h}/{n_t} transcripts) — complete knockout")
                elif lbl == "SELECTIVE":
                    pts.append(f"⚠ Isoform: SELECTIVE ({n_h}/{n_t}) — verify target isoform")
                else:
                    pts.append(f"✗ Isoform: SINGLE ({n_h}/{n_t}) — incomplete knockdown risk")

            gnom = bd.get("gnomad") or c.get("gnomad")
            if gnom and gnom.get("gnomad_checked"):
                risk = gnom.get("risk_level", "NONE")
                af   = gnom.get("af_max", 0) or 0
                if risk == "NONE":
                    pts.append("✅ No seed-region SNPs detected (gnomAD v4)")
                elif risk in ("LOW", "MODERATE"):
                    pts.append(f"⚠ Minor population variant overlap (gnomAD risk: {risk})")
                else:
                    pts.append(f"✗ Common seed-region SNP (AF {af:.3f}) — binding disrupted in diverse samples")

            if ot_sc >= 75:
                pts.append(f"✅ Low predicted off-target risk (score {ot_sc}/100)")
            elif ot_sc >= 50:
                pts.append(f"⚠ Moderate off-target risk ({ot_sc}/100) — await CRISPOR genome-wide results")
            else:
                pts.append(f"✗ Higher off-target risk ({ot_sc}/100) — validate with GUIDE-seq before use")

            if rank_i == 0:
                pts.append("🏆 Highest composite score across all evaluated criteria in this run")
            elif rank_i <= 2:
                pts.append(f"📊 Ranked #{rank_i+1} — strong alternative; include in experimental shortlist")

            c["explainability"] = pts

            # Do-not-use flags
            dnu_reasons: list[str] = []
            if gc > 0.75:
                dnu_reasons.append(f"GC {gc_pct}% > 75% — secondary structure risk")
            if gc < 0.35:
                dnu_reasons.append(f"GC {gc_pct}% < 35% — insufficient melting temp")
            for base in "ACGT":
                if base * 5 in seq:
                    dnu_reasons.append(f"Homopolymer run ({base}×5) — truncates Pol III transcription")
                    break
            if "TTTTT" in seq:
                dnu_reasons.append("Poly-T → Pol III termination signal")
            if gnom and gnom.get("risk_level") == "HIGH":
                dnu_reasons.append("Seed-region SNP >1% AF — Cas9 binding disrupted")
            if c.get("safety_label") == "AVOID":
                dnu_reasons.append("Safety tier AVOID — high off-target prediction")

            c["do_not_use"]  = len(dnu_reasons) > 0
            c["dnu_reasons"] = dnu_reasons

        # ── Rejected guide examples (transparency) ────────────────────────
        # The pipeline already ran; expose a note about rejection in output.
        # (Actual rejected objects not preserved in JSON path; use metadata.)
        meta = output.get("metadata", {})
        n_rejected = (meta.get("total_candidates_evaluated", 0) -
                      meta.get("candidates_passed_filters", 0))
        output["rejection_summary"] = {
            "n_rejected": n_rejected,
            "note": (
                f"{n_rejected} candidate guides were rejected before ranking. "
                "Common reasons: GC content outside 40–70%, homopolymer runs (≥5 same bases), "
                "PAM site quality, or seed-region GC > 75%. "
                "Rejection logic is intentional — unsafe or low-quality guides are excluded "
                "to protect result reliability."
            ),
        }

        # ── Audit trail ───────────────────────────────────────────────────
        run_ts = datetime.now(timezone.utc)
        output["audit_trail"] = {
            "run_id":    f"SNIP-{run_ts.strftime('%Y%m%d')}-{uuid.uuid4().hex[:6].upper()}",
            "timestamp": run_ts.isoformat(),
            "input_gene":     gene or "custom_sequence",
            "organism":       organism,
            "cas_variant":    cas_variant,
            "guide_length":   guide_length,
            "model_version":  "azimuth-gbr-v1.0 (GBR, n_estimators=400, trained on Doench 2016)",
            "database_versions": {
                "gnomad":           "v4 (GRCh38) — queried live",
                "clinvar":          "current (NCBI Entrez API — queried live)",
                "refseq":           "current (NCBI Entrez API — queried live)",
                "crispor_genome":   "hg38" if organism == "human" else organism,
                "azimuth_training": "Doench et al. 2016 (n=5,310 experimentally measured guides)",
            },
            "parameters": {
                "min_gc":  min_gc,
                "max_gc":  max_gc,
                "top_n":   top_n,
            },
            "known_limitations": [
                "gnomAD coordinate mapping uses linear mRNA→genomic scaling — may be inaccurate for highly intronic genes",
                "Off-target score before CRISPOR returns is a heuristic (seed GC + self-complementarity), not genome-wide",
                "mRNA positions reported, not GRCh38 genomic coordinates",
                "Model trained on SpCas9 data — off-target heuristics less calibrated for other Cas variants",
                "Isoform check uses string matching in transcript sequences — does not account for RNA folding",
            ],
        }

        return output

    finally:
        tmp_path.unlink(missing_ok=True)


def _run_batch_pipeline(
    gene_list: list[str],
    organism: str,
    cas_variant: str,
    top_n: int,
) -> dict:
    """
    Batch design: fetch sequences for each gene and run pipeline.
    Returns combined results with per-gene candidate lists.
    """
    results: dict[str, dict] = {}

    for gene in gene_list[:10]:   # hard cap at 10 genes per batch
        gene = gene.strip().upper()
        if not gene:
            continue
        try:
            fasta_text, accession = _fetch_sequence_entrez(gene, organism)
            fasta_bytes = fasta_text.encode()
            gene_result = _run_pipeline(
                fasta_bytes=fasta_bytes,
                cas_variant=cas_variant,
                guide_length=20,
                min_gc=0.40,
                max_gc=0.70,
                top_n=top_n,
                organism=organism,
                gene_symbol=gene,
            )
            gene_result["gene"] = gene
            gene_result["accession"] = accession
            results[gene] = gene_result
        except Exception as exc:
            results[gene] = {"gene": gene, "error": str(exc), "candidates": []}

    return {"batch": True, "genes": list(results.keys()), "results": results}


# ── Routes ─────────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def root():
    return (_static / "index.html").read_text()


@app.get("/variants")
async def list_variants():
    return {
        v: {"pattern": cfg["pattern"], "position": cfg["position"]}
        for v, cfg in PAM_REGISTRY.items()
    }


@app.get("/fetch-gene")
async def fetch_gene(
    gene: str = Query(...),
    organism: str = Query("human"),
):
    gene = gene.strip()
    if not gene or len(gene) > 50:
        raise HTTPException(400, "Invalid gene name")
    try:
        fasta_text, accession = _fetch_sequence_entrez(gene, organism)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    except Exception as exc:
        raise HTTPException(500, f"Unexpected error: {exc}")

    seq_len = len("".join(
        l.strip() for l in fasta_text.splitlines() if not l.startswith(">")
    ))
    return JSONResponse({
        "gene": gene, "organism": organism,
        "accession": accession, "fasta": fasta_text, "length": seq_len,
    })


@app.post("/design")
async def design(
    file: UploadFile = File(...),
    cas_variant: str = Query("SpCas9"),
    guide_length: int = Query(20, ge=17, le=25),
    min_gc: float = Query(0.40, ge=0.0, le=1.0),
    max_gc: float = Query(0.70, ge=0.0, le=1.0),
    top_n: int = Query(20, ge=1, le=200),
    organism: str = Query("human"),
    gene_symbol: str = Query("", description="Gene symbol for gnomAD/isoform/ClinVar annotation"),
):
    """
    Accept a FASTA upload, queue the pipeline, return a job_id immediately.
    Client polls GET /job/{job_id} for status + results.
    """
    if cas_variant not in PAM_REGISTRY:
        raise HTTPException(400, f"Unknown Cas variant '{cas_variant}'")
    if min_gc >= max_gc:
        raise HTTPException(400, "min_gc must be less than max_gc")

    fasta_bytes = await file.read()
    if len(fasta_bytes) == 0:
        raise HTTPException(400, "Empty file uploaded")

    text_preview = fasta_bytes[:200].decode("utf-8", errors="replace")
    if not any(c in text_preview for c in (">", "A", "C", "G", "T", "a", "c", "g", "t")):
        raise HTTPException(400, "File does not appear to be a valid FASTA")

    job_id = queue.submit(
        _run_pipeline,
        fasta_bytes, cas_variant, guide_length, min_gc, max_gc, top_n, organism, gene_symbol,
    )

    return JSONResponse({"job_id": job_id, "status": "queued"}, status_code=202)


@app.post("/batch-design")
async def batch_design(
    genes: str = Query(..., description="Comma-separated gene symbols (max 10)"),
    organism: str = Query("human"),
    cas_variant: str = Query("SpCas9"),
    top_n: int = Query(5, ge=1, le=20),
):
    """
    Batch design mode: design guides for multiple genes at once.
    Accepts comma-separated gene symbols, returns job_id.
    """
    gene_list = [g.strip() for g in genes.split(",") if g.strip()][:10]
    if not gene_list:
        raise HTTPException(400, "No valid gene symbols provided")
    if cas_variant not in PAM_REGISTRY:
        raise HTTPException(400, f"Unknown Cas variant '{cas_variant}'")

    job_id = queue.submit(
        _run_batch_pipeline,
        gene_list, organism, cas_variant, top_n,
    )

    return JSONResponse({
        "job_id": job_id,
        "status": "queued",
        "genes": gene_list,
        "mode": "batch",
    }, status_code=202)


@app.get("/job/{job_id}")
async def job_status(job_id: str):
    """
    Poll for pipeline job status.

    Response shapes:
      {"status": "queued",  "progress": "Queued…",  "elapsed_s": 0.1}
      {"status": "running", "progress": "Scoring…", "elapsed_s": 3.4}
      {"status": "done",    "result": {...},         "elapsed_s": 8.1}
      {"status": "failed",  "error": "...",          "elapsed_s": 2.0}
    """
    job = queue.get(job_id)
    if not job:
        raise HTTPException(404, f"Job '{job_id}' not found (may have expired)")
    return JSONResponse(job.to_dict())


@app.get("/crispor-scores")
async def crispor_scores(
    batch_id: str = Query(...),
    gene_symbol: str = Query("", description="Gene symbol for ClinVar off-target annotation"),
):
    """Poll CRISPOR for real off-target results, annotated with ClinVar data."""
    if not batch_id or not batch_id.replace("-", "").replace("_", "").isalnum() or len(batch_id) > 30:
        raise HTTPException(400, "Invalid batch_id")

    scores = crispor_fetch(batch_id)
    if scores is None:
        return JSONResponse({"status": "pending"})

    converted = {}
    for seq, data in scores.items():
        entry = {**data, "snipgen_offtarget_score": crispor_to_offtarget_score(data)}

        # ClinVar annotation: annotate the gene locus of off-target hits
        locus = data.get("gene_locus", "")
        if locus and locus != "—":
            # locus format: "exon:GENENAME" or "intron:GENENAME"
            parts = locus.split(":")
            gene_name = parts[-1].strip() if len(parts) >= 2 else ""
            if gene_name:
                try:
                    gene_ann = clinvar_annotate_gene(gene_name)
                    entry["clinvar_offtarget"] = {
                        "gene":   gene_name,
                        "tier":   gene_ann.get("tier", "MINIMAL"),
                        "variants": gene_ann.get("variants", 0),
                        "disease": gene_ann.get("disease", ""),
                        "color":  gene_ann.get("color", "#9ca3af"),
                        "label":  gene_ann.get("label", ""),
                    }
                except Exception:
                    pass

        converted[seq] = entry

    return JSONResponse({"status": "ready", "scores": converted})


@app.get("/cloning-primers")
async def cloning_primers(
    guide: str = Query(..., description="20-mer guide sequence (no PAM)"),
    vector: str = Query("pX330"),
):
    """Generate cloning oligos for a single guide + vector combination."""
    guide = guide.strip().upper()
    if len(guide) < 17 or len(guide) > 25:
        raise HTTPException(400, "Guide must be 17-25 nt")
    if any(n not in "ACGTN" for n in guide):
        raise HTTPException(400, "Guide must contain only ACGTN")
    return JSONResponse(design_cloning_oligos(guide, vector))


@app.get("/base-edit")
async def base_edit(
    guide: str = Query(..., description="20-mer guide sequence (no PAM)"),
    strand: str = Query("+"),
):
    """Analyze a guide for base editing suitability (CBE/ABE)."""
    guide = guide.strip().upper()
    if len(guide) < 17:
        raise HTTPException(400, "Guide must be ≥17 nt")
    return JSONResponse(analyze_base_editing(guide, strand))
