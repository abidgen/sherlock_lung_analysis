# Sherlock-Lung Presentation Analysis
**Abid Al Reza, PhD**  

## Scientific Question

How does the somatic mutational landscape differ between never-smoker and smoker lung adenocarcinoma, and what does this tell us about endogenous mutational processes, driver biology, immune context, and copy-number support in never-smoker lung cancer?

## Dataset

This project uses public TCGA LUAD data.

Required local input files:

```text
data/TCGA_LUAD/
├── TCGA_LUAD_somatic.maf.gz
├── TCGA_LUAD_HiSeqV2.gz
└── TCGA_LUAD_clinical.tsv
```

Optional local input file for integrated CNA analysis:

```text
data/TCGA_LUAD/
└── TCGA_PANCAN_CNA.tsv.gz      # or TCGA_LUAD_CNA.tsv / TCGA_LUAD_CNA.tsv.gz
```

Data sources:

- Somatic mutation calls: TCGA LUAD MAF from GDC/API-based download.
- Expression matrix: UCSC Xena TCGA LUAD HiSeqV2.
- Clinical annotation: UCSC Xena TCGA LUAD phenotype/clinical matrix.
- Optional CNA: TCGA PanCancer or LUAD GISTIC2 thresholded copy-number matrix.

## Relevance to Sherlock-Lung

The Sherlock-Lung studies characterized the somatic landscape of never-smoker lung cancer at large scale. This analysis independently reproduces and extends related observations using public TCGA LUAD data:

1. Lower tobacco-associated mutational signal in never-smokers.
2. EGFR-dominant driver profile in never-smoker LUAD.
3. Distinct tumor microenvironment and immune-score patterns.
4. Differential expression and Hallmark pathway activation by smoking group.
5. Integrated EGFR mutation, expression, and optional copy-number support.

## Project Structure

```text
sherlock_lung/
├── scripts/
│   ├── 01_TCGA_LUAD_analysis.R       # Expression, ESTIMATE, Hallmark/KEGG GSEA, survival
│   ├── 02_MAF_analysis.R             # Somatic mutation / MAF analysis
│   ├── 03_integrated_analysis.R      # EGFR expression, EGFR-mutant DE, optional CNA
│   └── run_analysis.sh               # Full pipeline wrapper
├── data/TCGA_LUAD/                   # Input data; not committed to git
├── figures/                          # Generated figures; not committed to git
├── logs/                             # Timestamped run logs; not committed to git
├── envs/
│   └── sherlock_lung_environment.yml # Conda environment pin
├── README.md
└── .gitignore
```

## Reproducibility

This project uses a **conda-only reproducibility setup**. The conda environment file pins the R runtime plus CRAN and Bioconductor dependencies, including `maftools`, `limma`, `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`, `msigdbr`, `estimate`, `ggplot2`, `dplyr`, `ggpubr`, `survminer`, and related packages.

Primary reproducibility file:

```text
envs/sherlock_lung_environment.yml
```

Minimum files to keep under version control:

```text
envs/sherlock_lung_environment.yml
scripts/01_TCGA_LUAD_analysis.R
scripts/02_MAF_analysis.R
scripts/03_integrated_analysis.R
scripts/run_analysis.sh
README.md
.gitignore
```

Generated outputs and large inputs should stay out of git:

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

`renv` is intentionally not used because it conflicted with conda-installed Bioconductor packages by overriding `.libPaths()`.

## Environment Setup

Create the conda environment from the project root:

```bash
conda env create -f envs/sherlock_lung_environment.yml
conda activate sherlock_lung
```

If the environment already exists:

```bash
conda activate sherlock_lung
```

To update the environment pin after package changes:

```bash
mkdir -p envs
conda env export --no-builds > envs/sherlock_lung_environment.yml
```

To confirm key packages are available:

```bash
Rscript -e 'pkgs <- c("maftools", "limma", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "msigdbr", "estimate", "ggplot2", "dplyr", "ggpubr", "survminer"); sapply(pkgs, function(p) { library(p, character.only=TRUE); as.character(packageVersion(p)) })'
```

## Running the Analysis

Run the full pipeline from the project root:

```bash
conda activate sherlock_lung
chmod +x scripts/run_analysis.sh
./scripts/run_analysis.sh
```

The wrapper executes:

1. `scripts/02_MAF_analysis.R` — required.
2. `scripts/01_TCGA_LUAD_analysis.R` — required.
3. `scripts/03_integrated_analysis.R` — optional; the wrapper continues if this step fails due to missing optional CNA data.

The run wrapper also checks that the three required input files are present before starting the analysis.

Timestamped logs are written to:

```text
logs/maf_run_YYYYMMDD_HHMMSS.log
logs/expression_run_YYYYMMDD_HHMMSS.log
logs/integrated_run_YYYYMMDD_HHMMSS.log
logs/full_run_YYYYMMDD_HHMMSS.log
```

Run individual scripts manually if needed:

```bash
Rscript scripts/02_MAF_analysis.R
Rscript scripts/01_TCGA_LUAD_analysis.R
Rscript scripts/03_integrated_analysis.R
```

## Analysis Scripts

- `scripts/02_MAF_analysis.R` performs somatic mutation analysis, including MAF summaries, oncoplot, Ti/Tv comparison, mutational signatures, COSMIC_4/SBS4 contribution comparison, driver comparison, EGFR lollipop plot, EGFR mutation label table, optional maftools pathway summaries, and a summary statistics table.
- `scripts/01_TCGA_LUAD_analysis.R` performs limma differential expression, ESTIMATE immune/stromal scoring, MSigDB Hallmark GSEA when `msigdbr` is available, KEGG fallback otherwise, and immune-score survival analysis.
- `scripts/03_integrated_analysis.R` performs integrated driver-expression analysis, EGFR-mutant vs wildtype expression analysis in never-smokers, and optional CNA analysis using `TCGA_LUAD_CNA.tsv`, `TCGA_LUAD_CNA.tsv.gz`, or `TCGA_PANCAN_CNA.tsv.gz`.

## Figures and Tables Generated

### Expression / immune / pathway / survival

| Figure | Output file | Description | Sherlock-Lung connection |
|---|---|---|---|
| Fig1 | `Fig1_Volcano_NeverSmoker_vs_Smoker.png` | Differential expression volcano plot | Expression differences by smoking status |
| Fig2 | `Fig2_ESTIMATE_Immune_Scores.png` | ESTIMATE immune, stromal, and composite scores | Tumor microenvironment differences |
| Fig3 | `Fig3_GSEA_Hallmark.png` or `Fig3_GSEA_KEGG.png` | GSEA dotplot | Differential pathway activation |
| Fig4 | `Fig4_Survival_ImmuneScore.png` | Overall survival by immune-score group, if enough events are available | Immune context and prognosis |

### MAF / mutation-landscape analysis

| Figure/table | Output file | Description | Sherlock-Lung connection |
|---|---|---|---|
| MAF Fig1 | `MAF_Fig1_Summary_NeverSmoker.png` | Never-smoker mutation summary | Lower / distinct mutational burden |
| MAF Fig1b | `MAF_Fig1b_Summary_Smoker.png` | Smoker mutation summary | Smoker comparison baseline |
| MAF Fig2 | `MAF_Fig2_Oncoplot_NeverSmoker.png` | Top mutated genes in never-smokers | EGFR-dominant never-smoker profile |
| MAF Fig3 | `MAF_Fig3_TiTv_Comparison.png` | Ti/Tv comparison | Smoking-related mutational process differences |
| MAF Fig4 | `MAF_Fig4_Signatures_NeverSmoker.png` | Never-smoker mutational signatures | Low tobacco-signature contribution |
| MAF Fig4b | `MAF_Fig4b_Signatures_Smoker.png` | Smoker mutational signatures | Tobacco-associated comparison |
| MAF Fig4c | `MAF_Fig4c_Contributions_NeverSmoker.png` | Per-sample signature contribution | Within-group signature heterogeneity |
| MAF Fig4d | `MAF_Fig4d_COSMIC4_Contribution_Comparison.png` | COSMIC_4/SBS4 contribution comparison | Tobacco signature lower in never-smokers |
| MAF Fig5 | `MAF_Fig5_Driver_Comparison.png` | Driver-gene forest plot | EGFR vs smoking-associated drivers |
| MAF Fig6 | `MAF_Fig6_EGFR_Lollipop.png` | EGFR lollipop plot | EGFR hotspot characterization |
| MAF Fig6b | `MAF_Fig6b_EGFR_Label_Table.csv` / `.png` | EGFR mutation label table | Slide-friendly mutation labels |
| MAF Fig7 | `MAF_Fig7_OncogenicPathways_NeverSmoker.png` | Optional maftools oncogenic pathway plot, if supported by installed maftools | Pathway-level driver summary |
| MAF Fig7b | `MAF_Fig7b_OncogenicPathways_Smoker.png` | Optional maftools oncogenic pathway plot for smokers | Smoker pathway comparison |
| Summary | `Summary_Statistics_Table.png` and `summary_statistics.csv` | Slide-friendly numeric summary | Key mutation statistics |

### Integrated EGFR / CNA analysis

| Figure/table | Output file | Description | Requirement |
|---|---|---|---|
| Fig5 | `Fig5_EGFR_Expression_by_Mutation.png` | EGFR expression by mutation status and smoking group | Required core files |
| Fig5b | `Fig5b_EGFR_Mutant_vs_Wildtype_Volcano.png` | EGFR-mutant vs wildtype DE in never-smokers, if sample counts are sufficient | Required core files |
| Fig7 | `Fig7_CNA_Landscape.png` | Amplification/deletion frequencies for selected LUAD genes | Optional CNA file |
| Fig7b | `Fig7b_EGFR_Mutation_Amplification.png` | EGFR mutation plus amplification status by smoking group | Optional CNA file |
| CNA table | `CNA_summary.csv` | Gene-level CNA summary | Optional CNA file |

## Notes on Current Implementation

- The expression analysis uses **limma**, not DESeq2, because the UCSC Xena HiSeqV2 matrix is already normalized/log-transformed rather than raw counts.
- ESTIMATE uses the GCT writer that is compatible with the installed `estimate` package. ESTIMATE sample IDs are normalized before clinical matching to avoid TCGA barcode issues caused by R name conversion.
- GSEA uses MSigDB Hallmark gene sets through `msigdbr` when available. The code uses the current `msigdbr` `ncbi_gene` column and falls back to KEGG GSEA if `msigdbr` is unavailable.
- Survival event coding uses TCGA-style clinical values such as `DECEASED` and `LIVING`; survival figures are skipped if too few events are available.
- The MAF pathway section uses `maftools::pathways()` / `maftools::plotPathways()` only if those functions exist in the installed maftools version.
- The integrated CNA analysis is optional and supports either LUAD-specific CNA files or the PanCancer file `TCGA_PANCAN_CNA.tsv.gz`.
- The expression matrix is filtered to primary tumor samples only (TCGA barcode positions 14–15 = "01") before differential expression, ESTIMATE, GSEA, and survival analysis, removing normal adjacent tissue (type 11) and recurrent tumor (type 02) samples. The filtered dataset includes 515 primary tumor samples (75 never-smoker, 422 smoker, remainder unannotated for smoking status).
- Telomere analysis was removed from the main workflow because the LUAD-specific Xena telomere path was not reliably downloadable.

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


