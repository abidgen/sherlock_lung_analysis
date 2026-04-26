# Sherlock-Lung Presentation Analysis
**Abid Al Reza, PhD**  
Prepared for: Computational Scientist II interview, CGR/FNL/NCI (req4512)

## Scientific Question
How does the somatic mutational landscape differ between never-smoker and 
smoker lung adenocarcinoma, and what does this tell us about endogenous 
mutational processes in never-smoker lung cancer?

## Dataset
- **TCGA LUAD** (The Cancer Genome Atlas, Lung Adenocarcinoma)
- 440 samples with somatic WGS/WES mutation calls
- Clinical annotation including tobacco smoking history
- Source: GDC API (somatic MAF) + UCSC Xena (expression + clinical)

## Relevance to Sherlock-Lung
The Sherlock-Lung cohort (Zhang et al. Nature Genetics 2021; Díaz-Gay et al. 
Nature 2025) characterized the somatic landscape of never-smoker lung cancer 
at unprecedented scale. This analysis independently reproduces and extends 
key findings from those papers using public TCGA data:

1. Absence of SBS4 (tobacco) signature in never-smokers
2. EGFR dominance in never-smoker LUAD
3. Distinct tumor microenvironment immune landscape
4. Differential pathway activation between groups

## Reproducibility
- All analysis containerized in Docker (rocker/r-ver:4.3.1)
- R package versions pinned in Dockerfile
- Data downloaded from public repositories (GDC, UCSC Xena)
- All scripts version-controlled

## Environment
- Docker: rocker/r-ver:4.3.1
- R: 4.3.1
- Bioconductor: 3.17
- Key packages: maftools 2.16.0, Seurat 4.3.0, harmony 0.1.0

## Analysis Scripts
- `scripts/02_MAF_analysis.R` — Somatic mutation landscape
- `scripts/01_TCGA_LUAD_analysis.R` — Differential expression + immune

## Figures Generated
| Figure | Description | Sherlock-Lung Connection |
|--------|-------------|--------------------------|
| MAF_Fig1 | Mutation summary never-smoker vs smoker | Low TMB in never-smokers |
| MAF_Fig2 | Oncoplot top 20 genes never-smoker | EGFR dominance |
| MAF_Fig3 | TiTv comparison | Distinct mutational processes |
| MAF_Fig4 | Trinucleotide spectrum | SBS4 absent in never-smokers |
| MAF_Fig5 | Driver gene forest plot | EGFR vs KRAS enrichment |
| MAF_Fig6 | EGFR lollipop plot | Hotspot characterization |

## Data Sources
- Somatic MAF: https://api.gdc.cancer.gov/data/c06465a3-50e7-46f7-b2dd-7bd654ca206b
- Expression: UCSC Xena TCGA LUAD HiSeqV2
- Clinical: UCSC Xena TCGA LUAD clinical matrix
