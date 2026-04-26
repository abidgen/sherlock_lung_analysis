# Sherlock-Lung Presentation Analysis

**Abid Al Reza, PhD**  
Prepared for: Computational Scientist II interview, CGR/FNL/NCI (req4512)

## Scientific Question

How does the somatic mutational landscape differ between never-smoker and smoker lung adenocarcinoma, and what does this tell us about endogenous mutational processes in never-smoker lung cancer?

## Dataset

- **TCGA LUAD** (The Cancer Genome Atlas, Lung Adenocarcinoma)
- Somatic mutation calls from TCGA LUAD MAF
- Clinical annotation including tobacco smoking history
- Expression matrix from UCSC Xena TCGA LUAD HiSeqV2
- Source: GDC API for somatic MAF + UCSC Xena for expression and clinical data

Expected local data layout:

```text
data/TCGA_LUAD/
├── TCGA_LUAD_somatic.maf.gz
├── TCGA_LUAD_HiSeqV2.gz
└── TCGA_LUAD_clinical.tsv
```

## Relevance to Sherlock-Lung

The Sherlock-Lung studies characterized the somatic landscape of never-smoker lung cancer at large scale. This analysis independently reproduces and extends related observations using public TCGA LUAD data:

1. Lower tobacco-associated mutational signal in never-smokers
2. EGFR dominance in never-smoker LUAD
3. Distinct tumor microenvironment / immune-score patterns
4. Differential pathway activation between never-smoker and smoker LUAD

## Project Structure

```text
sherlock_lung/
├── scripts/
│   ├── 01_TCGA_LUAD_analysis.R      # Expression, ESTIMATE, GSEA, survival
│   ├── 02_MAF_analysis.R            # Somatic mutation / MAF analysis
│   └── run_analysis.sh              # Full reproducible run wrapper
├── data/TCGA_LUAD/                  # Input data; not committed to git
├── figures/                         # Generated figures; not committed to git
├── logs/                            # Timestamped run logs; not committed to git
├── envs/
│   └── sherlock_lung_environment.yml
├── README.md
└── .gitignore
```

## Reproducibility

This project uses a **conda-only reproducibility setup**. The environment file pins the R runtime plus the CRAN and Bioconductor packages required for the analysis, including `maftools`, `limma`, `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`, `ggplot2`, `dplyr`, `ggpubr`, `survminer`, and related dependencies.

The primary reproducibility file is:

```text
envs/sherlock_lung_environment.yml
```

Minimum files to keep under version control:

```text
envs/sherlock_lung_environment.yml
scripts/01_TCGA_LUAD_analysis.R
scripts/02_MAF_analysis.R
scripts/run_analysis.sh
README.md
.gitignore
```

Generated outputs and large input files should stay out of git:

```text
data/
figures/
logs/
```

Recommended `.gitignore`:

```text
data/
figures/
logs/
renv/
renv.lock
.Rprofile
.Rhistory
.RData
.Ruserdata
```

## Environment Setup

From the project root:

```bash
conda env create -f envs/sherlock_lung_environment.yml
conda activate sherlock_lung
```

If the conda environment already exists:

```bash
conda activate sherlock_lung
```

To update the environment pin after changing packages:

```bash
mkdir -p envs
conda env export --no-builds > envs/sherlock_lung_environment.yml
```

To confirm key packages are available:

```bash
Rscript -e 'pkgs <- c("maftools", "limma", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "ggplot2", "dplyr", "ggpubr", "survminer"); sapply(pkgs, function(p) { library(p, character.only=TRUE); as.character(packageVersion(p)) })'
```

## Running the Analysis

Run the full pipeline from the project root:

```bash
conda activate sherlock_lung
chmod +x scripts/run_analysis.sh
./scripts/run_analysis.sh
```

The run wrapper executes:

1. `scripts/02_MAF_analysis.R`
2. `scripts/01_TCGA_LUAD_analysis.R`

It writes timestamped logs to:

```text
logs/maf_run_YYYYMMDD_HHMMSS.log
logs/expression_run_YYYYMMDD_HHMMSS.log
logs/full_run_YYYYMMDD_HHMMSS.log
```

Run individual scripts manually if needed:

```bash
Rscript scripts/02_MAF_analysis.R
Rscript scripts/01_TCGA_LUAD_analysis.R
```

## Analysis Scripts

- `scripts/02_MAF_analysis.R` — somatic mutation landscape, MAF summaries, oncoplot, Ti/Tv comparison, mutational signatures, driver comparison, and EGFR lollipop plot.
- `scripts/01_TCGA_LUAD_analysis.R` — differential expression with limma, ESTIMATE immune/stromal scores, KEGG GSEA, and immune-score survival analysis.

## Figures Generated

| Figure | Output file | Description | Sherlock-Lung connection |
|---|---|---|---|
| Fig1 | `Fig1_Volcano_NeverSmoker_vs_Smoker.png` | Differential expression volcano plot | Expression differences by smoking status |
| Fig2 | `Fig2_ESTIMATE_Immune_Scores.png` | ESTIMATE immune, stromal, and composite scores | Tumor microenvironment differences |
| Fig3 | `Fig3_GSEA_Pathways.png` | KEGG GSEA dotplot | Differential pathway activation |
| Fig4 | `Fig4_Survival_ImmuneScore.png` | Overall survival by immune-score group | Immune context and prognosis |
| MAF Fig1 | `MAF_Fig1_Summary_NeverSmoker.png` | Never-smoker mutation summary | Lower / distinct mutational burden |
| MAF Fig1b | `MAF_Fig1b_Summary_Smoker.png` | Smoker mutation summary | Smoker comparison baseline |
| MAF Fig2 | `MAF_Fig2_Oncoplot_NeverSmoker.png` | Top mutated genes in never-smokers | EGFR-dominant never-smoker profile |
| MAF Fig3 | `MAF_Fig3_TiTv_Comparison.png` | Ti/Tv comparison | Smoking-related mutational process differences |
| MAF Fig4 | `MAF_Fig4_Signatures_NeverSmoker.png` | Never-smoker mutational signatures | Low tobacco-signature contribution |
| MAF Fig4b | `MAF_Fig4b_Signatures_Smoker.png` | Smoker mutational signatures | Tobacco-associated comparison |
| MAF Fig4c | `MAF_Fig4c_Contributions_NeverSmoker.png` | Per-sample signature contribution | Within-group signature heterogeneity |
| MAF Fig4d | `MAF_Fig4d_COSMIC4_Contribution_Comparison.png` | COSMIC_4/SBS4 contribution comparison | Tobacco signature lower in never-smokers |
| MAF Fig5 | `MAF_Fig5_Driver_Comparison.png` | Filtered driver-gene forest plot | EGFR vs smoking-associated drivers |
| MAF Fig6 | `MAF_Fig6_EGFR_Lollipop.png` | EGFR lollipop plot | EGFR hotspot characterization |
| MAF Fig6b | `MAF_Fig6b_EGFR_Label_Table.csv/.png` | EGFR mutation label table | Slide-friendly mutation labels |

## Notes on Current Implementation

- The expression analysis uses **limma**, not DESeq2, because the UCSC Xena HiSeqV2 matrix is already normalized/log-transformed rather than raw counts.
- ESTIMATE sample IDs are normalized before clinical matching to avoid TCGA barcode issues caused by R name conversion.
- Survival analysis uses `days_to_death` when available and falls back to `days_to_last_followup` for censored samples.
- The MAF driver comparison figure is filtered to significant genes with enough mutated samples to keep the forest plot readable.
- The EGFR lollipop plot labels only key interpretable positions; the full EGFR mutation label table is exported separately for slides.
- `renv` is intentionally not used because it conflicted with conda-installed Bioconductor packages by overriding `.libPaths()`.

## Data Sources

- Somatic MAF: GDC API TCGA LUAD somatic mutation file
- Expression: UCSC Xena TCGA LUAD HiSeqV2
- Clinical: UCSC Xena TCGA LUAD clinical matrix

## Recreating a Run

A minimal reproducible run should include:

```bash
conda env create -f envs/sherlock_lung_environment.yml
conda activate sherlock_lung
./scripts/run_analysis.sh
```

After completion, check:

```bash
ls -lh figures/*.png
ls -lh logs/*.log
```
