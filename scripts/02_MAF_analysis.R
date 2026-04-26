# ============================================================
# TCGA LUAD Somatic MAF Analysis
# Never-Smoker vs Smoker — Mutational Landscape
# Abid Al Reza, PhD — Sherlock-Lung Presentation
# ============================================================

BASE    <- "/analysis"
DATA    <- file.path(BASE, "data/TCGA_LUAD")
FIGURES <- file.path(BASE, "figures")
dir.create(FIGURES, showWarnings=FALSE)

library(maftools)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ---- 1. LOAD CLINICAL AND DEFINE SMOKING GROUPS ----
message("Loading clinical data...")
clin <- read.delim(file.path(DATA, "TCGA_LUAD_clinical.tsv"),
                   header=TRUE, sep="\t")

# tobacco_smoking_history: 1=never, 2=current, 3-4=former
clin_smoking <- clin %>%
  mutate(
    patient_id = substr(sampleID, 1, 12),
    smoking_group = case_when(
      tobacco_smoking_history == 1 ~ "Never_Smoker",
      tobacco_smoking_history %in% c(2,3,4) ~ "Smoker",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(smoking_group)) %>%
  select(patient_id, smoking_group, tobacco_smoking_history)

message(paste("Never-smokers:", sum(clin_smoking$smoking_group=="Never_Smoker")))
message(paste("Smokers:", sum(clin_smoking$smoking_group=="Smoker")))

# ---- 2. LOAD MAF ----
message("Loading MAF file...")
luad_maf <- read.maf(
  maf = file.path(DATA, "TCGA_LUAD_somatic.maf.gz"),
  clinicalData = clin_smoking %>% 
    rename(Tumor_Sample_Barcode=patient_id),
  verbose = FALSE
)

message("MAF loaded successfully.")
luad_maf

# ---- 3. SUBSET BY SMOKING STATUS ----
never_barcodes <- clin_smoking %>%
  filter(smoking_group=="Never_Smoker") %>%
  pull(patient_id)

smoker_barcodes <- clin_smoking %>%
  filter(smoking_group=="Smoker") %>%
  pull(patient_id)

# Match to full barcodes in MAF
maf_barcodes <- getSampleSummary(luad_maf)$Tumor_Sample_Barcode
never_maf_bc <- maf_barcodes[substr(maf_barcodes,1,12) %in% never_barcodes]
smoker_maf_bc <- maf_barcodes[substr(maf_barcodes,1,12) %in% smoker_barcodes]

message(paste("Never-smoker samples in MAF:", length(never_maf_bc)))
message(paste("Smoker samples in MAF:", length(smoker_maf_bc)))

maf_never  <- subsetMaf(luad_maf, tsb=never_maf_bc)
maf_smoker <- subsetMaf(luad_maf, tsb=smoker_maf_bc)

# ---- FIGURE 1 — MUTATION SUMMARY ----
message("Generating Figure 1: MAF Summary...")
png(file.path(FIGURES, "MAF_Fig1_Summary_NeverSmoker.png"),
    width=1400, height=700, res=150)
plotmafSummary(maf_never,
               rmOutlier=TRUE,
               addStat="median",
               dashboard=TRUE,
               titvRaw=FALSE,
               top=10,
               fs=0.8)
title("Never-Smoker LUAD — Somatic Mutation Summary (TCGA)",
      line=2, cex.main=1)
dev.off()

png(file.path(FIGURES, "MAF_Fig1b_Summary_Smoker.png"),
    width=1400, height=700, res=150)
plotmafSummary(maf_smoker,
               rmOutlier=TRUE,
               addStat="median",
               dashboard=TRUE,
               titvRaw=FALSE,
               top=10,
               fs=0.8)
title("Smoker LUAD — Somatic Mutation Summary (TCGA)",
      line=2, cex.main=1)
dev.off()
message("Figure 1 saved.")

# ---- FIGURE 2 — ONCOPLOT NEVER SMOKERS ----
message("Generating Figure 2: Oncoplot...")
png(file.path(FIGURES, "MAF_Fig2_Oncoplot_NeverSmoker.png"),
    width=1400, height=900, res=150)
oncoplot(maf_never,
         top=20,
         fontSize=10,
         legendFontSize=10,
         showTumorSampleBarcodes=FALSE,
         drawRowBar=TRUE,
         drawColBar=TRUE,
         title="Never-Smoker LUAD — Top 20 Mutated Genes (TCGA, n=never_smoker_n)")
dev.off()
message("Figure 2 saved.")

# ---- FIGURE 3 — MUTATION SPECTRUM COMPARISON ----
message("Generating Figure 3: Mutation spectrum...")
png(file.path(FIGURES, "MAF_Fig3_TiTv_Comparison.png"),
    width=1200, height=600, res=150)
par(mfrow=c(1,2))

titv_never  <- titv(maf=maf_never,  plot=FALSE, useSyn=TRUE)
titv_smoker <- titv(maf=maf_smoker, plot=FALSE, useSyn=TRUE)

plotTiTv(res=titv_never,
         colorCode=c("C>A"="#E41A1C",
                     "C>G"="#377EB8",
                     "C>T"="#4DAF4A",
                     "T>A"="#984EA3",
                     "T>C"="#FF7F00",
                     "T>G"="#A65628"),
         showBarcodes=FALSE)
title("Never-Smoker", cex.main=1.2)

plotTiTv(res=titv_smoker,
         colorCode=c("C>A"="#E41A1C",
                     "C>G"="#377EB8",
                     "C>T"="#4DAF4A",
                     "T>A"="#984EA3",
                     "T>C"="#FF7F00",
                     "T>G"="#A65628"),
         showBarcodes=FALSE)
title("Smoker", cex.main=1.2)

mtext("Transition/Transversion Comparison — Never-Smoker vs Smoker LUAD",
      outer=TRUE, cex=1.2, line=-1)
dev.off()
message("Figure 3 saved.")

# ---- FIGURE 4 — TRINUCLEOTIDE CONTEXT ----
message("Generating Figure 4: Trinucleotide mutation spectrum...")

# This is the key figure — shows absence of SBS4 (C>A tobacco) in never-smokers
png(file.path(FIGURES, "MAF_Fig4_Trinucleotide_NeverSmoker.png"),
    width=1400, height=600, res=150)
plotSignatures(maf_never,
               contribution=TRUE,
               title_size=1,
               show_title=TRUE)
title("Never-Smoker LUAD — Mutational Signatures (TCGA)",
      line=2.5, cex.main=1)
dev.off()

png(file.path(FIGURES, "MAF_Fig4b_Trinucleotide_Smoker.png"),
    width=1400, height=600, res=150)
plotSignatures(maf_smoker,
               contribution=TRUE,
               title_size=1,
               show_title=TRUE)
title("Smoker LUAD — Mutational Signatures (TCGA)",
      line=2.5, cex.main=1)
dev.off()
message("Figure 4 saved.")

# ---- FIGURE 5 — DRIVER GENE COMPARISON ----
message("Generating Figure 5: Driver gene comparison...")

# Compare mutation frequencies between groups
comp <- mafCompare(m1=maf_never,
                   m2=maf_smoker,
                   m1Name="Never_Smoker",
                   m2Name="Smoker",
                   minMut=5)

png(file.path(FIGURES, "MAF_Fig5_Driver_Comparison.png"),
    width=1000, height=800, res=150)
forestPlot(mafCompareRes=comp,
           pVal=0.05,
           geneFontSize=0.8,
           titleSize=0.8)
dev.off()
message("Figure 5 saved.")

# ---- FIGURE 6 — LOLLIPOP EGFR ----
message("Generating Figure 6: EGFR lollipop...")
png(file.path(FIGURES, "MAF_Fig6_EGFR_Lollipop.png"),
    width=1200, height=600, res=150)
lollipopPlot(maf=maf_never,
             gene="EGFR",
             AACol="HGVSp_Short",
             showMutationRate=TRUE,
             labelPos="all",
             title="EGFR Mutations — Never-Smoker LUAD (TCGA)")
dev.off()
message("Figure 6 saved.")

# ---- SUMMARY ----
message("============================================")
message("MAF analysis complete. Figures generated:")
list.files(FIGURES, pattern="MAF_")
message("============================================")

# Key statistics for presentation
message("\n=== KEY STATISTICS FOR PRESENTATION ===")
message(paste("Never-smoker median mutations:", 
              median(getSampleSummary(maf_never)$total)))
message(paste("Smoker median mutations:", 
              median(getSampleSummary(maf_smoker)$total)))
message(paste("Top gene in never-smokers:", 
              getGeneSummary(maf_never)$Hugo_Symbol[1],
              "(",
              round(getGeneSummary(maf_never)$MutatedSamples[1]/
                    length(never_maf_bc)*100, 1),
              "%)"))
