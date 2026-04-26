# ============================================================
# TCGA LUAD Somatic MAF Analysis
# Never-Smoker vs Smoker — Mutational Landscape
# Abid Al Reza, PhD — Sherlock-Lung Presentation
# ============================================================

BASE    <- "/media/wrath/bioinfor_learning/sherlock_lung"
DATA    <- file.path(BASE, "data/TCGA_LUAD")
FIGURES <- file.path(BASE, "figures")
dir.create(FIGURES, showWarnings=FALSE)

library(maftools)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ---- 1. LOAD CLINICAL AND DEFINE SMOKING GROUPS ----
# Xena clinical sampleIDs: TCGA-05-4244-01 (15 chars)
# MAF barcodes:            TCGA-05-4244-01A-11D-A24D-08 (full)
# Match key: first 12 chars (TCGA-TSS-Patient, e.g. TCGA-05-4244)
message("Loading clinical data...")
clin <- read.delim(file.path(DATA, "TCGA_LUAD_clinical.tsv"),
                   header=TRUE, sep="\t")

# tobacco_smoking_history: 1=never, 2=current, 3-4=former
clin_smoking <- clin %>%
  mutate(
    patient_id    = substr(sampleID, 1, 12),   # 12-char patient key
    smoking_group = case_when(
      tobacco_smoking_history == 1             ~ "Never_Smoker",
      tobacco_smoking_history %in% c(2, 3, 4) ~ "Smoker",
      TRUE                                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(smoking_group)) %>%
  select(patient_id, smoking_group, tobacco_smoking_history)

message(paste("Never-smokers:", sum(clin_smoking$smoking_group == "Never_Smoker")))
message(paste("Smokers:",       sum(clin_smoking$smoking_group == "Smoker")))

# ---- 2. LOAD MAF ----
message("Loading MAF file...")
luad_maf_raw <- read.maf(
  maf     = file.path(DATA, "TCGA_LUAD_somatic.maf.gz"),
  verbose = FALSE
)
message("MAF loaded successfully.")

# ---- 3. SUBSET BY SMOKING STATUS ----
# Build a barcode → patient_id map and join with clinical smoking groups
maf_barcodes <- as.character(getSampleSummary(luad_maf_raw)$Tumor_Sample_Barcode)

# One aliquot per patient — keep first (sorted lexicographically, so 01A < 01B)
maf_annotated <- data.frame(
  Tumor_Sample_Barcode = maf_barcodes,
  patient_id           = substr(maf_barcodes, 1, 12),
  stringsAsFactors     = FALSE
) %>%
  arrange(Tumor_Sample_Barcode) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  inner_join(clin_smoking, by = "patient_id")

never_maf_bc  <- maf_annotated %>% filter(smoking_group == "Never_Smoker") %>% pull(Tumor_Sample_Barcode)
smoker_maf_bc <- maf_annotated %>% filter(smoking_group == "Smoker")       %>% pull(Tumor_Sample_Barcode)

message(paste("Never-smoker samples in MAF:", length(never_maf_bc)))
message(paste("Smoker samples in MAF:",       length(smoker_maf_bc)))

# Attach clinical annotations to MAF object (needed for maftools plots)
clin_for_maf <- maf_annotated %>%
  select(Tumor_Sample_Barcode, smoking_group, tobacco_smoking_history)

luad_maf <- read.maf(
  maf          = file.path(DATA, "TCGA_LUAD_somatic.maf.gz"),
  clinicalData = clin_for_maf,
  verbose      = FALSE
)

maf_never  <- subsetMaf(luad_maf, tsb=never_maf_bc)
maf_smoker <- subsetMaf(luad_maf, tsb=smoker_maf_bc)

# ---- FIGURE 1 — MUTATION SUMMARY ----
message("Generating Figure 1: MAF Summary...")

png(file.path(FIGURES, "MAF_Fig1_Summary_NeverSmoker.png"),
    width=2400, height=1150, res=200)

par(oma=c(0, 0, 3, 0))  # add outer top margin

plotmafSummary(
  maf_never,
  rmOutlier=TRUE,
  addStat="median",
  dashboard=TRUE,
  titvRaw=FALSE,
  top=10,
  fs=1.0,
  titleSize=0.5
)

mtext("Never-Smoker LUAD — Somatic Mutation Summary (TCGA)",
      side=3, outer=TRUE, line=1, cex=1.1, font=2)

dev.off()


png(file.path(FIGURES, "MAF_Fig1b_Summary_Smoker.png"),
    width=2000, height=1150, res=200)

par(oma=c(0, 0, 3, 0))  # add outer top margin

plotmafSummary(
  maf_smoker,
  rmOutlier=TRUE,
  addStat="median",
  dashboard=TRUE,
  titvRaw=FALSE,
  top=10,
  fs=1.0,
  titleSize=0.5
)

mtext("Smoker LUAD — Somatic Mutation Summary (TCGA)",
      side=3, outer=TRUE, line=1, cex=1.1, font=2)

dev.off()

message("Figure 1 saved.")

# ---- FIGURE 2 — ONCOPLOT NEVER SMOKERS ----
message("Generating Figure 2: Oncoplot...")
png(file.path(FIGURES, "MAF_Fig2_Oncoplot_NeverSmoker.png"),
    width=1600, height=1000, res=200)
oncoplot(maf_never,
         top=20,
         fontSize=0.9,
         legendFontSize=1.1,
         titleFontSize=1.2,
         showTumorSampleBarcodes=FALSE,
         drawRowBar=TRUE,
         drawColBar=TRUE,
         titleText=paste0("Never-Smoker LUAD — Top 20 Mutated Genes (TCGA, n=", length(never_maf_bc), ")"))
dev.off()
message("Figure 2 saved.")

# ---- FIGURE 3 — MUTATION SPECTRUM COMPARISON ----
# plotTiTv calls graphics::layout() internally, which overrides par(mfrow).
# Fix: render each group to a temp PNG, then tile them with grid::rasterGrob.
message("Generating Figure 3: Mutation spectrum...")

library(grid)
library(gridExtra)

titv_col <- c("C>A"="#E41A1C", "C>G"="#377EB8", "C>T"="#4DAF4A",
              "T>A"="#984EA3", "T>C"="#FF7F00", "T>G"="#A65628")

titv_never  <- titv(maf=maf_never,  plot=FALSE, useSyn=TRUE)
titv_smoker <- titv(maf=maf_smoker, plot=FALSE, useSyn=TRUE)

# Render each panel to its own temp file at the final resolution
tmp_never  <- tempfile(fileext=".png")
tmp_smoker <- tempfile(fileext=".png")

png(tmp_never, width=900, height=700, res=200)
plotTiTv(res=titv_never, color=titv_col, showBarcodes=FALSE)
title("Never-Smoker", cex.main=1.4, font.main=2)
dev.off()

png(tmp_smoker, width=900, height=700, res=200)
plotTiTv(res=titv_smoker, color=titv_col, showBarcodes=FALSE)
title("Smoker", cex.main=1.4, font.main=2)
dev.off()

img_never  <- png::readPNG(tmp_never)
img_smoker <- png::readPNG(tmp_smoker)

png(file.path(FIGURES, "MAF_Fig3_TiTv_Comparison.png"),
    width=2000, height=800, res=200)
grid.arrange(
  rasterGrob(img_never,  interpolate=TRUE),
  rasterGrob(img_smoker, interpolate=TRUE),
  ncol=2,
  top=textGrob("Transition/Transversion Comparison — Never-Smoker vs Smoker LUAD",
               gp=gpar(fontsize=13, fontface="bold"))
)
dev.off()
unlink(c(tmp_never, tmp_smoker))
message("Figure 3 saved.")

# ---- FIGURE 4 — MUTATIONAL SIGNATURES ----
# maftools 2.x workflow: trinucleotideMatrix → extractSignatures → plotSignatures
# Key figure: shows absence of SBS4 (C>A tobacco) in never-smokers
message("Generating Figure 4: Mutational signatures...")

sig_rds_never  <- file.path(DATA, "sig_never.rds")
sig_rds_smoker <- file.path(DATA, "sig_smoker.rds")

library(BSgenome.Hsapiens.UCSC.hg38)
library(NMF)  # must load after BSgenome to restore NMF's seed generic

if (!file.exists(sig_rds_never)) {
  mat_never <- trinucleotideMatrix(maf_never,  ref_genome="BSgenome.Hsapiens.UCSC.hg38",
                                   prefix=NULL, add=TRUE, useSyn=TRUE)
  sig_never <- extractSignatures(mat_never, n=6, parallel=4)
  saveRDS(sig_never, sig_rds_never)
} else {
  sig_never <- readRDS(sig_rds_never)
}

if (!file.exists(sig_rds_smoker)) {
  mat_smoker <- trinucleotideMatrix(maf_smoker, ref_genome="BSgenome.Hsapiens.UCSC.hg38",
                                    prefix=NULL, add=TRUE, useSyn=TRUE)
  sig_smoker <- extractSignatures(mat_smoker, n=6, parallel=4)
  saveRDS(sig_smoker, sig_rds_smoker)
} else {
  sig_smoker <- readRDS(sig_rds_smoker)
}

png(file.path(FIGURES, "MAF_Fig4_Signatures_NeverSmoker.png"),
    width=1800, height=1000, res=200)
plotSignatures(nmfRes=sig_never,
               contributions=FALSE,
               title_size=1.0,
               font_size=0.7,
               show_title=TRUE)
dev.off()

png(file.path(FIGURES, "MAF_Fig4b_Signatures_Smoker.png"),
    width=1800, height=1000, res=200)
plotSignatures(nmfRes=sig_smoker,
               contributions=FALSE,
               title_size=1.0,
               font_size=0.7,
               show_title=TRUE)
dev.off()
message("Figure 4 saved.")

# ---- FIGURE 4c — PER-SAMPLE CONTRIBUTION (NEVER-SMOKER) ----
# plotSignatures(contributions=TRUE) calls graphics::layout() internally, so
# it needs a standalone PNG with enough width for 69 sample columns.
message("Generating Figure 4c: Per-sample contribution plot (never-smoker)...")
png(file.path(FIGURES, "MAF_Fig4c_Contributions_NeverSmoker.png"),
    width=3400, height=900, res=200)
plotSignatures(nmfRes=sig_never,
               contributions=TRUE,
               title_size=0.85,
               font_size=0.55,
               show_title=TRUE,
               show_barcodes=FALSE)
dev.off()
message("Figure 4c saved.")

# ---- FIGURE 4d — COSMIC_4 CONTRIBUTION: NEVER-SMOKER vs SMOKER ----
# Cosine similarity to COSMIC_4 is high for one never-smoker signature (0.958)
# but the actual per-sample contribution fraction is low relative to smokers.
message("Generating Figure 4d: COSMIC_4 contribution comparison...")

cs_n <- compareSignatures(sig_never,  sig_db="legacy")$cosine_similarities
cs_s <- compareSignatures(sig_smoker, sig_db="legacy")$cosine_similarities

# Identify extracted signatures whose best COSMIC match is COSMIC_4
best_n <- apply(cs_n, 1, function(x) names(which.max(x)))
best_s <- apply(cs_s, 1, function(x) names(which.max(x)))
c4_sigs_n <- names(best_n[best_n == "COSMIC_4"])  # e.g. Signature_4
c4_sigs_s <- names(best_s[best_s == "COSMIC_4"])  # e.g. Signature_3, Signature_4

message(paste("Never-smoker COSMIC_4-matched:", paste(c4_sigs_n, collapse=", "),
              "| cosine:", paste(round(cs_n[c4_sigs_n, "COSMIC_4"], 3), collapse="/")))
message(paste("Smoker COSMIC_4-matched:", paste(c4_sigs_s, collapse=", "),
              "| cosine:", paste(round(cs_s[c4_sigs_s, "COSMIC_4"], 3), collapse="/")))

# Sum contributions across all COSMIC_4-matched signatures within each group
contrib_n_c4 <- colSums(sig_never$contributions[c4_sigs_n,  , drop=FALSE])
contrib_s_c4 <- colSums(sig_smoker$contributions[c4_sigs_s, , drop=FALSE])

cosine_n_str <- paste(round(cs_n[c4_sigs_n, "COSMIC_4"], 3), collapse="/")
cosine_s_str <- paste(round(cs_s[c4_sigs_s, "COSMIC_4"], 3), collapse="/")

never_label  <- paste0("Never-Smoker\n(n=", length(contrib_n_c4), ", cosine=", cosine_n_str, ")")
smoker_label <- paste0("Smoker\n(n=", length(contrib_s_c4), ", cosine=", cosine_s_str, ")")

dat4d <- rbind(
  data.frame(group=never_label,  frac_pct=contrib_n_c4 * 100),
  data.frame(group=smoker_label, frac_pct=contrib_s_c4 * 100)
)
dat4d$group <- factor(dat4d$group, levels=c(never_label, smoker_label))

pval <- wilcox.test(contrib_n_c4, contrib_s_c4)$p.value
plab <- ifelse(pval < 0.001, "Wilcoxon p < 0.001", sprintf("Wilcoxon p = %.3f", pval))
y_annot <- max(dat4d$frac_pct) * 1.07

p4d <- ggplot(dat4d, aes(x=group, y=frac_pct, fill=group)) +
  geom_violin(trim=FALSE, alpha=0.50, color=NA) +
  geom_boxplot(width=0.10, outlier.shape=16, outlier.size=1.2,
               fill="white", color="grey30", linewidth=0.5) +
  geom_jitter(width=0.07, alpha=0.22, size=0.9, color="grey30") +
  scale_fill_manual(values=c("#2166AC", "#D73027")) +
  annotate("segment", x=1, xend=2,
           y=y_annot * 1.01, yend=y_annot * 1.01, color="grey50") +
  annotate("text", x=1.5, y=y_annot * 1.07, label=plab, size=4.5) +
  labs(
    title="COSMIC_4 (Tobacco / SBS4) Signature: Contribution vs Cosine Similarity",
    subtitle=paste0(
      "Cosine similarity to COSMIC_4 is high (never ", cosine_n_str,
      " vs smoker ", cosine_s_str, "), yet the actual per-sample\n",
      "contribution fraction is significantly lower in never-smokers — ",
      "reflecting minimal tobacco mutagenesis"
    ),
    x=NULL,
    y="% of mutations attributed to COSMIC_4-matched signature(s)"
  ) +
  theme_classic(base_size=13) +
  theme(legend.position="none",
        plot.title    = element_text(face="bold", size=14),
        plot.subtitle = element_text(size=10.5, color="grey35", lineheight=1.3),
        axis.text.x   = element_text(size=12),
        axis.text.y   = element_text(size=11))

ggsave(file.path(FIGURES, "MAF_Fig4d_COSMIC4_Contribution_Comparison.png"),
       p4d, width=8, height=6.5, dpi=200)
message("Figure 4d saved.")

# ---- FIGURE 5 — DRIVER GENE COMPARISON ----
message("Generating Figure 5: Driver gene comparison...")

comp <- mafCompare(m1=maf_never,
                   m2=maf_smoker,
                   m1Name="Never_Smoker",
                   m2Name="Smoker",
                   minMut=3)

# Show top 20 by p-value regardless of significance threshold
comp_top <- comp
comp_top$results <- comp$results[order(comp$results$pval)][1:min(20, nrow(comp$results))]

png(file.path(FIGURES, "MAF_Fig5_Driver_Comparison.png"),
    width=2000, height=1000, res=200)
forestPlot(mafCompareRes=comp_top,
           pVal=1.0,
           geneFontSize=1.0,
           titleSize=1.1)
dev.off()
message("Figure 5 saved.")

# ---- FIGURE 6 — LOLLIPOP EGFR ----
message("Generating Figure 6: EGFR lollipop...")

label_positions <- c(62, 746, 750, 858, 861)

png(file.path(FIGURES, "MAF_Fig6_EGFR_Lollipop.png"),
    width=2400, height=1200, res=200)

par(oma=c(0, 0, 3, 0))
par(mar=c(7, 5, 4, 2))
par(las=1)

lollipopPlot(
  maf=maf_never,
  gene="EGFR",
  AACol="HGVSp_Short",
  showMutationRate=TRUE,
  labelPos=label_positions,
  repel=TRUE,
  domainLabelSize=0.01,
  axisTextSize=c(1.0, 1.0),
  legendTxtSize=0.9,
  titleSize=c(1.2, 0.9)
)

mtext("EGFR Mutations — Never-Smoker LUAD (TCGA)",
      side=3, outer=TRUE, line=1, cex=1.1, font=2)

dev.off()

message("Figure 6 saved.")

# ---- FIGURE 6b — EGFR MUTATION LABEL TABLE ----
message("Generating EGFR mutation label table...")

egfr_label_table <- maf_never@data %>%
  dplyr::filter(Hugo_Symbol == "EGFR", !is.na(HGVSp_Short)) %>%
  dplyr::mutate(
    aa_pos = as.numeric(gsub(".*?(\\d+).*", "\\1", HGVSp_Short))
  ) %>%
  dplyr::filter(!is.na(aa_pos)) %>%
  dplyr::group_by(aa_pos, HGVSp_Short, Variant_Classification) %>%
  dplyr::summarise(
    Mutated_Samples = dplyr::n_distinct(Tumor_Sample_Barcode),
    .groups = "drop"
  ) %>%
  dplyr::arrange(aa_pos, HGVSp_Short) %>%
  dplyr::rename(
    Position = aa_pos,
    Label = HGVSp_Short,
    Mutation_Type = Variant_Classification
  )

# Save as CSV for editing in PowerPoint/Excel
write.csv(
  egfr_label_table,
  file.path(FIGURES, "MAF_Fig6b_EGFR_Label_Table.csv"),
  row.names = FALSE
)

# Save as PNG table for slide
library(gridExtra)
library(grid)

png(file.path(FIGURES, "MAF_Fig6b_EGFR_Label_Table.png"),
    width = 1800, height = 1200, res = 200)

grid.newpage()

table_plot <- tableGrob(
  egfr_label_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    core = list(fg_params = list(hjust = 0, x = 0.05)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

grid.draw(table_plot)

dev.off()

message("EGFR mutation label table saved.")


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
