# SnipGen

**AI-driven CRISPR guide RNA design platform**

[![Live Demo](https://img.shields.io/badge/live-snipgen--1.onrender.com-brightgreen)](https://snipgen-1.onrender.com)
[![Python](https://img.shields.io/badge/python-3.10+-blue)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

SnipGen takes a gene name, fetches the mRNA sequence from NCBI automatically, finds all valid guide RNAs, ranks them using a machine learning model trained on real experimental data, and returns everything a researcher needs to order and clone those guides — all in one browser session, no login required.

---

## Live Demo

[snipgen-1.onrender.com](https://snipgen-1.onrender.com)

---

## Features

### Core Pipeline

- **Gene auto-fetch**: type a gene symbol (e.g. TP53, BRCA1, KRAS) and SnipGen pulls the official RefSeq mRNA from NCBI automatically. Supports human, mouse, zebrafish, rat, Drosophila, and C. elegans. No FASTA upload required.
- **Multi-Cas variant support**: design guides for SpCas9 (NGG), SaCas9 (NNGRRT), Cpf1/Cas12a (TTTN), xCas9, and Cas9-NG. PAM scanning and filters adjust per variant.
- **GC content filtering, homopolymer removal, and deduplication** applied to all candidates before scoring.

### ML On-Target Scoring

SnipGen uses a GradientBoostingRegressor trained on **5,310 real CRISPR screen guides** from the Doench 2016 / Azimuth dataset — the same experimental data used to train the published Azimuth model (Nature Biotechnology).

- Feature vector: 97 dimensions (nucleotide identity × position, GC content, melting temperature, homopolymer runs, seed-region GC, di-nucleotide features)
- Spearman correlation on held-out validation: **0.556** (published Azimuth: ~0.56–0.60)
- Model file: `snipgen/scoring/models/ontarget_xgb.pkl`, retrained on every deployment

### Genome-Wide Off-Target Analysis (CRISPOR)

After returning initial results, SnipGen submits your sequence to [CRISPOR](http://crispor.gi.ucsc.edu) (Haussler lab, UCSC), which runs Cas-OFFinder against the full reference genome (hg38/mm39). Each guide card updates live with:

- MIT specificity score (0–100)
- CFD off-target score
- Mismatch site counts (0/1/2/3 mismatches)

> Note: the initial off-target score shown before CRISPOR returns is a seed-region heuristic (GC content + self-complementarity), not a genome search result. It is labeled as an estimate.

### ClinVar Off-Target Consequence Annotation

When CRISPOR identifies an off-target site in a gene, SnipGen cross-references that gene against a curated database of ~150 clinically significant genes (COSMIC tier 1/2 cancer genes + ACMG SF v3.2 actionable genes). Off-target hits in high-consequence genes (e.g. BRCA2, PTEN) are flagged CRITICAL in red. No other free CRISPR tool does this.

### Cloning Primer Auto-Design

For every guide, SnipGen generates the annealing oligos for direct cloning:

```
Top oligo:    CACC + [G if needed] + guide_20mer
Bottom oligo: AAAC + reverse_complement(guide) + [C if G added]
```

Supported vectors: pX330, pX458, pX459 (BbsI), lentiCRISPRv2 (BsmBI), pGuide (Esp3I), pX601/AAV (BsaI), T7 in vitro transcription. Click any oligo to copy to clipboard. Includes Tm display and ordering notes.

### Base Editing Compatibility

Analyzes each guide for CBE and ABE suitability. Highlights the edit window (positions 4–8) per-nucleotide: C in window = CBE targetable, A in window = ABE targetable. Detects bystander edits. Efficiency estimate based on Rees & Liu 2018 position data.

### Isoform Specificity Analysis

Fetches all RefSeq transcript variants for the target gene and checks each guide against every transcript. Labels guides as: PAN-ISOFORM (all transcripts), BROAD (>75%), SELECTIVE (25–75%), or SINGLE-ISOFORM.

### gnomAD Population Variant Filter

For human guides, queries gnomAD v4 for common population variants (>1% allele frequency) overlapping each guide's seed region. Flags guides where a common SNP could disrupt Cas9 binding.

> Known limitation: genomic coordinate mapping uses linear scaling from mRNA position, which does not account for introns. gnomAD lookups may be offset for multi-exon genes. Fix in progress.

### Batch Gene Design

Enter up to 10 gene symbols at once. SnipGen runs the full pipeline for each in parallel and returns a combined results view.

### Guide Comparison Modal

Select up to 4 guides to compare side-by-side: sequence, GC%, on-target score, off-target score, safety label, isoform coverage, ClinVar flags, and cloning oligos. Best value per row is highlighted.

---

## Tech Stack

| Layer | Stack |
|---|---|
| Backend | Python, Flask, BioPython |
| ML | scikit-learn (GradientBoostingRegressor) |
| Off-target | CRISPOR (external API), Cas-OFFinder |
| Annotations | NCBI Entrez, gnomAD v4 API, ClinVar (curated local DB) |
| Frontend | Vanilla JS, HTML/CSS |
| Deployment | Render (persistent server, async job queue) |
| CI | GitHub Actions |

---

## Installation

```bash
pip install biopython
pip install -e ".[dev]"
```

## Usage

```bash
# Design gRNAs from a FASTA file
snipgen design --input target.fasta --output-dir results/

# With custom options
snipgen design --input target.fasta \
               --output-dir results/ \
               --format csv json \
               --cas-variant SpCas9 \
               --guide-length 20 \
               --min-gc 0.40 \
               --max-gc 0.70 \
               --top-n 20 \
               --verbose

# Validate input only
snipgen validate --input target.fasta

# List supported Cas variants
snipgen list-variants
```

## ML Model Integration

Drop in a trained sklearn model (joblib-serialized):

```bash
snipgen design --input target.fasta --ml-model model.joblib --ml-weight 0.4
```

The model receives an 84-dimensional feature vector per candidate (20-pos one-hot + 4 scalar features).

---

## Tests

```bash
pytest tests/
```

Test coverage includes: CLI interface, FASTA reader, GC filter, ML scorer, off-target filter, PAM filter, full pipeline, and rule scorer.

---

## Known Limitations

- gnomAD coordinate mapping uses linear mRNA-to-genomic scaling; does not account for intron positions. Affects variant lookups on multi-exon genes.
- Off-target score displayed before CRISPOR returns is a heuristic proxy, not genome-wide data.
- Output coordinates are mRNA positions (RefSeq accession), not GRCh38 chr:pos.

---

## Project Background

Built as an independent biomedical engineering project to explore practical CRISPR tool development, ML-based guide scoring, and full-stack bioinformatics deployment. The ML model is trained from scratch on the Doench 2016 Azimuth dataset at deployment time, not a wrapper around existing tools.

---

## Author

**Lakshya Dharwal**
B.S.E. Biomedical Engineering, Arizona State University (May 2026)
[github.com/lakshya-dharwal](https://github.com/lakshya-dharwal) | ldharwal2003@gmail.com
