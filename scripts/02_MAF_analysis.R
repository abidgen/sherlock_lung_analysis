# ============================================================
# TCGA LUAD Somatic MAF Analysis
# Never-Smoker vs Smoker — Mutational Landscape
# Abid Al Reza, PhD — Sherlock-Lung Presentation
# ============================================================

# ---- Paths ----
BASE <- Sys.getenv("SHERLOCK_LUNG_BASE", unset = getwd())
DATA <- file.path(BASE, "data", "TCGA_LUAD")
FIGURES <- file.path(BASE, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----
required_pkgs <- c("maftools", "dplyr", "ggplot2", "grid", "gridExtra", "png")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them in the conda environment and rerun.")
}

library(maftools)
library(dplyr)
library(ggplot2)

# ---- Helpers ----
check_file <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) stop("Missing or empty file: ", path)
  invisible(path)
}
patient12 <- function(x) substr(gsub("\\.", "-", x), 1, 12)
safe_pct <- function(x, denom) if (length(x) == 0 || denom == 0) NA_real_ else as.numeric(x) / denom * 100

# ---- Input files ----
clin_file <- check_file(file.path(DATA, "TCGA_LUAD_clinical.tsv"))
maf_file <- check_file(file.path(DATA, "TCGA_LUAD_somatic.maf.gz"))

# ============================================================
# 1. Load clinical data and define smoking groups
# ============================================================
message("Loading clinical data...")
clin <- read.delim(clin_file, header = TRUE, sep = "\t")
if (!"sampleID" %in% colnames(clin)) colnames(clin)[1] <- "sampleID"

clin_smoking <- clin %>%
  mutate(
    patient_id = patient12(sampleID),
    smoking_group = case_when(
      tobacco_smoking_history == 1 ~ "Never_Smoker",
      tobacco_smoking_history %in% c(2, 3, 4) ~ "Smoker",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(smoking_group)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  select(patient_id, smoking_group, tobacco_smoking_history)

message("Never-smokers in clinical: ", sum(clin_smoking$smoking_group == "Never_Smoker"))
message("Smokers in clinical: ", sum(clin_smoking$smoking_group == "Smoker"))

# ============================================================
# 2. Load MAF and subset by smoking group
# ============================================================
message("Loading MAF file...")
luad_maf_raw <- read.maf(maf = maf_file, verbose = FALSE)
message("MAF loaded successfully.")

maf_barcodes <- as.character(getSampleSummary(luad_maf_raw)$Tumor_Sample_Barcode)
maf_annotated <- data.frame(
  Tumor_Sample_Barcode = maf_barcodes,
  patient_id = patient12(maf_barcodes),
  stringsAsFactors = FALSE
) %>%
  arrange(Tumor_Sample_Barcode) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  inner_join(clin_smoking, by = "patient_id")

never_maf_bc <- maf_annotated %>% filter(smoking_group == "Never_Smoker") %>% pull(Tumor_Sample_Barcode)
smoker_maf_bc <- maf_annotated %>% filter(smoking_group == "Smoker") %>% pull(Tumor_Sample_Barcode)

message("Never-smoker samples in MAF: ", length(never_maf_bc))
message("Smoker samples in MAF: ", length(smoker_maf_bc))
if (length(never_maf_bc) < 2 || length(smoker_maf_bc) < 2) stop("Too few matched MAF samples.")

clin_for_maf <- maf_annotated %>% select(Tumor_Sample_Barcode, smoking_group, tobacco_smoking_history)
luad_maf <- read.maf(maf = maf_file, clinicalData = clin_for_maf, verbose = FALSE)
maf_never <- subsetMaf(luad_maf, tsb = never_maf_bc)
maf_smoker <- subsetMaf(luad_maf, tsb = smoker_maf_bc)

# ============================================================
# 3. Mutation summary plots
# ============================================================
message("Generating Figure 1: MAF summary plots...")

png(file.path(FIGURES, "MAF_Fig1_Summary_NeverSmoker.png"), width = 2400, height = 1150, res = 200)
par(oma = c(0, 0, 3, 0))
plotmafSummary(maf_never, rmOutlier = TRUE, addStat = "median", dashboard = TRUE,
               titvRaw = FALSE, top = 10, fs = 1.0, titleSize = 0.5)
mtext("Never-Smoker LUAD — Somatic Mutation Summary (TCGA)", side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2)
dev.off()

png(file.path(FIGURES, "MAF_Fig1b_Summary_Smoker.png"), width = 2000, height = 1150, res = 200)
par(oma = c(0, 0, 3, 0))
plotmafSummary(maf_smoker, rmOutlier = TRUE, addStat = "median", dashboard = TRUE,
               titvRaw = FALSE, top = 10, fs = 1.0, titleSize = 0.5)
mtext("Smoker LUAD — Somatic Mutation Summary (TCGA)", side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2)
dev.off()
message("Figure 1 saved.")

# ============================================================
# 4. Oncoplot and Ti/Tv comparison
# ============================================================
message("Generating Figure 2: Oncoplot...")
png(file.path(FIGURES, "MAF_Fig2_Oncoplot_NeverSmoker.png"), width = 1600, height = 1000, res = 200)
oncoplot(maf_never, top = 20, fontSize = 0.9, legendFontSize = 1.1,
         titleFontSize = 1.2, showTumorSampleBarcodes = FALSE,
         drawRowBar = TRUE, drawColBar = TRUE,
         titleText = paste0("Never-Smoker LUAD — Top 20 Mutated Genes (TCGA, n=", length(never_maf_bc), ")"))
dev.off()
message("Figure 2 saved.")

message("Generating Figure 3: Ti/Tv comparison...")
library(grid)
library(gridExtra)

titv_col <- c("C>A" = "#E41A1C", "C>G" = "#377EB8", "C>T" = "#4DAF4A",
              "T>A" = "#984EA3", "T>C" = "#FF7F00", "T>G" = "#A65628")
titv_never <- titv(maf = maf_never, plot = FALSE, useSyn = TRUE)
titv_smoker <- titv(maf = maf_smoker, plot = FALSE, useSyn = TRUE)

tmp_never <- tempfile(fileext = ".png")
tmp_smoker <- tempfile(fileext = ".png")

png(tmp_never, width = 900, height = 700, res = 200)
plotTiTv(res = titv_never, color = titv_col, showBarcodes = FALSE)
title("Never-Smoker", cex.main = 1.4, font.main = 2)
dev.off()

png(tmp_smoker, width = 900, height = 700, res = 200)
plotTiTv(res = titv_smoker, color = titv_col, showBarcodes = FALSE)
title("Smoker", cex.main = 1.4, font.main = 2)
dev.off()

img_never <- png::readPNG(tmp_never)
img_smoker <- png::readPNG(tmp_smoker)
png(file.path(FIGURES, "MAF_Fig3_TiTv_Comparison.png"), width = 2000, height = 800, res = 200)
grid.arrange(rasterGrob(img_never, interpolate = TRUE),
             rasterGrob(img_smoker, interpolate = TRUE),
             ncol = 2,
             top = textGrob("Transition/Transversion Comparison — Never-Smoker vs Smoker LUAD",
                            gp = gpar(fontsize = 13, fontface = "bold")))
dev.off()
unlink(c(tmp_never, tmp_smoker))
message("Figure 3 saved.")

# ============================================================
# 5. Mutational signatures and COSMIC_4/SBS4 contribution
# ============================================================
message("Generating Figure 4: Mutational signatures...")
sig_rds_never <- file.path(DATA, "sig_never.rds")
sig_rds_smoker <- file.path(DATA, "sig_smoker.rds")

signature_available <- requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) &&
  requireNamespace("NMF", quietly = TRUE)

sig_never <- NULL
sig_smoker <- NULL
c4_mean_never <- NA_real_
c4_mean_smoker <- NA_real_

if (signature_available) {
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(NMF)

  if (!file.exists(sig_rds_never)) {
    mat_never <- trinucleotideMatrix(maf_never, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                                     prefix = NULL, add = TRUE, useSyn = TRUE)
    sig_never <- extractSignatures(mat_never, n = 6, parallel = 4)
    saveRDS(sig_never, sig_rds_never)
  } else {
    sig_never <- readRDS(sig_rds_never)
  }

  if (!file.exists(sig_rds_smoker)) {
    mat_smoker <- trinucleotideMatrix(maf_smoker, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                                      prefix = NULL, add = TRUE, useSyn = TRUE)
    sig_smoker <- extractSignatures(mat_smoker, n = 6, parallel = 4)
    saveRDS(sig_smoker, sig_rds_smoker)
  } else {
    sig_smoker <- readRDS(sig_rds_smoker)
  }

  png(file.path(FIGURES, "MAF_Fig4_Signatures_NeverSmoker.png"), width = 1800, height = 1000, res = 200)
  plotSignatures(nmfRes = sig_never, contributions = FALSE, title_size = 1.0, font_size = 0.7, show_title = TRUE)
  dev.off()

  png(file.path(FIGURES, "MAF_Fig4b_Signatures_Smoker.png"), width = 1800, height = 1000, res = 200)
  plotSignatures(nmfRes = sig_smoker, contributions = FALSE, title_size = 1.0, font_size = 0.7, show_title = TRUE)
  dev.off()

  png(file.path(FIGURES, "MAF_Fig4c_Contributions_NeverSmoker.png"), width = 3400, height = 900, res = 200)
  plotSignatures(nmfRes = sig_never, contributions = TRUE, title_size = 0.85,
                 font_size = 0.55, show_title = TRUE, show_barcodes = FALSE)
  dev.off()

  message("Generating Figure 4d: COSMIC_4 contribution comparison...")
  cs_n <- compareSignatures(sig_never, sig_db = "legacy")$cosine_similarities
  cs_s <- compareSignatures(sig_smoker, sig_db = "legacy")$cosine_similarities
  best_n <- apply(cs_n, 1, function(x) names(which.max(x)))
  best_s <- apply(cs_s, 1, function(x) names(which.max(x)))
  c4_sigs_n <- names(best_n[best_n == "COSMIC_4"])
  c4_sigs_s <- names(best_s[best_s == "COSMIC_4"])

  if (length(c4_sigs_n) > 0 && length(c4_sigs_s) > 0) {
    contrib_n_c4 <- colSums(sig_never$contributions[c4_sigs_n, , drop = FALSE])
    contrib_s_c4 <- colSums(sig_smoker$contributions[c4_sigs_s, , drop = FALSE])
    c4_mean_never <- mean(contrib_n_c4, na.rm = TRUE) * 100
    c4_mean_smoker <- mean(contrib_s_c4, na.rm = TRUE) * 100
    cosine_n_str <- paste(round(cs_n[c4_sigs_n, "COSMIC_4"], 3), collapse = "/")
    cosine_s_str <- paste(round(cs_s[c4_sigs_s, "COSMIC_4"], 3), collapse = "/")

    never_label <- paste0("Never-Smoker\n(n=", length(contrib_n_c4), ", cosine=", cosine_n_str, ")")
    smoker_label <- paste0("Smoker\n(n=", length(contrib_s_c4), ", cosine=", cosine_s_str, ")")
    dat4d <- rbind(data.frame(group = never_label, frac_pct = contrib_n_c4 * 100),
                   data.frame(group = smoker_label, frac_pct = contrib_s_c4 * 100))
    dat4d$group <- factor(dat4d$group, levels = c(never_label, smoker_label))
    pval <- wilcox.test(contrib_n_c4, contrib_s_c4)$p.value
    plab <- ifelse(pval < 0.001, "Wilcoxon p < 0.001", sprintf("Wilcoxon p = %.3f", pval))
    y_annot <- max(dat4d$frac_pct, na.rm = TRUE) * 1.07

    p4d <- ggplot(dat4d, aes(x = group, y = frac_pct, fill = group)) +
      geom_violin(trim = FALSE, alpha = 0.50, color = NA) +
      geom_boxplot(width = 0.10, outlier.shape = 16, outlier.size = 1.2,
                   fill = "white", color = "grey30", linewidth = 0.5) +
      geom_jitter(width = 0.07, alpha = 0.22, size = 0.9, color = "grey30") +
      scale_fill_manual(values = c("#2166AC", "#D73027")) +
      annotate("segment", x = 1, xend = 2, y = y_annot * 1.01, yend = y_annot * 1.01, color = "grey50") +
      annotate("text", x = 1.5, y = y_annot * 1.07, label = plab, size = 4.5) +
      labs(title = "COSMIC_4 (Tobacco / SBS4) Signature Contribution",
           subtitle = paste0("COSMIC_4 matched signatures: never cosine=", cosine_n_str,
                             " | smoker cosine=", cosine_s_str),
           x = NULL, y = "% of mutations attributed to COSMIC_4-matched signature(s)") +
      theme_classic(base_size = 13) +
      theme(legend.position = "none", plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 10.5, color = "grey35"))

    ggsave(file.path(FIGURES, "MAF_Fig4d_COSMIC4_Contribution_Comparison.png"),
           p4d, width = 8, height = 6.5, dpi = 200)
  } else {
    message("No COSMIC_4-matched signatures found in one or both groups; skipping Fig4d.")
  }
  message("Figure 4 saved.")
} else {
  message("Signature packages not available; skipping signature figures.")
}

# ============================================================
# 6. Driver comparison, EGFR lollipop, oncogenic pathways
# ============================================================
message("Generating Figure 5: Driver gene comparison...")
comp <- mafCompare(m1 = maf_never, m2 = maf_smoker,
                   m1Name = "Never_Smoker", m2Name = "Smoker", minMut = 5)
comp_plot <- comp
comp_plot$results <- comp$results[
  comp$results$pval < 0.05 & (comp$results$n1 >= 5 | comp$results$n2 >= 5),
]
if (nrow(comp_plot$results) == 0) {
  message("No significant genes after filtering; using top 20 by p-value for Fig5.")
  comp_plot$results <- comp$results[order(comp$results$pval), ][1:min(20, nrow(comp$results)), ]
  pval_cutoff <- 1.0
} else {
  comp_plot$results <- comp_plot$results[order(comp_plot$results$pval), ][1:min(25, nrow(comp_plot$results)), ]
  pval_cutoff <- 0.05
}

png(file.path(FIGURES, "MAF_Fig5_Driver_Comparison.png"), width = 2000, height = 1000, res = 200)
forestPlot(mafCompareRes = comp_plot, pVal = pval_cutoff, geneFontSize = 1.0, titleSize = 1.1)
dev.off()
message("Figure 5 saved.")

message("Generating Figure 6: EGFR lollipop...")
label_positions <- c(62, 746, 750, 858, 861)
png(file.path(FIGURES, "MAF_Fig6_EGFR_Lollipop.png"), width = 2400, height = 1200, res = 200)
par(oma = c(0, 0, 3, 0), mar = c(7, 5, 4, 2), las = 1)
lollipopPlot(maf = maf_never, gene = "EGFR", AACol = "HGVSp_Short",
             showMutationRate = TRUE, labelPos = label_positions, repel = TRUE,
             domainLabelSize = 0.01, axisTextSize = c(1.0, 1.0),
             legendTxtSize = 0.9, titleSize = c(1.2, 0.9))
mtext("EGFR Mutations — Never-Smoker LUAD (TCGA)", side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2)
dev.off()
message("Figure 6 saved.")

# ---- FIGURE 6b — EGFR MUTATION LABEL TABLE ----
message("Generating EGFR mutation label table...")

egfr_label_table <- maf_never@data %>%
  dplyr::filter(Hugo_Symbol == "EGFR", !is.na(HGVSp_Short), HGVSp_Short != "") %>%
  dplyr::mutate(
    Position = suppressWarnings(as.numeric(gsub(".*?(\\d+).*", "\\1", HGVSp_Short)))
  ) %>%
  dplyr::filter(!is.na(Position)) %>%
  dplyr::group_by(Position, HGVSp_Short, Variant_Classification) %>%
  dplyr::summarise(
    Mutated_Samples = dplyr::n_distinct(Tumor_Sample_Barcode),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Position, HGVSp_Short) %>%
  dplyr::rename(
    Label = HGVSp_Short,
    Mutation_Type = Variant_Classification
  )

write.csv(
  egfr_label_table,
  file.path(FIGURES, "MAF_Fig6b_EGFR_Label_Table.csv"),
  row.names = FALSE
)

library(gridExtra)
library(grid)

png(file.path(FIGURES, "MAF_Fig6b_EGFR_Label_Table.png"),
    width = 1800, height = 1200, res = 200)

grid.newpage()

table_plot <- gridExtra::tableGrob(
  egfr_label_table,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 9,
    core = list(fg_params = list(hjust = 0, x = 0.05)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

grid::grid.draw(table_plot)
dev.off()

message("EGFR mutation label table saved.")

# ---- FIGURE 7 — ONCOGENIC SIGNALING PATHWAYS ----
message("Generating Figure 7: Oncogenic signaling pathways...")

if (exists("pathways", where=asNamespace("maftools"), inherits=FALSE)) {

  pws_never <- maftools::pathways(
    maf = maf_never,
    plotType = "treemap"
  )

  pws_smoker <- maftools::pathways(
    maf = maf_smoker,
    plotType = "treemap"
  )

  # Save pathway summary tables if returned
  if (!is.null(pws_never)) {
    write.csv(
      pws_never,
      file.path(DATA, "oncogenic_pathways_never_smoker.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(pws_smoker)) {
    write.csv(
      pws_smoker,
      file.path(DATA, "oncogenic_pathways_smoker.csv"),
      row.names = FALSE
    )
  }

  # Plot pathway matrix, if available
  if (exists("plotPathways", where=asNamespace("maftools"), inherits=FALSE)) {
    png(file.path(FIGURES, "MAF_Fig7_OncogenicPathways_NeverSmoker.png"),
        width=1800, height=1200, res=200)
    maftools::plotPathways(maf = maf_never, pathlist = pws_never)
    dev.off()

    png(file.path(FIGURES, "MAF_Fig7b_OncogenicPathways_Smoker.png"),
        width=1800, height=1200, res=200)
    maftools::plotPathways(maf = maf_smoker, pathlist = pws_smoker)
    dev.off()

    message("Figure 7 pathway plots saved.")
  } else {
    message("plotPathways() not available; pathway tables saved only.")
  }

} else {
  message("Skipping Figure 7: maftools::pathways() not available in this maftools version.")
}

# ============================================================
# 7. Summary statistics table
# ============================================================
message("Generating summary statistics table...")
never_tmb <- getSampleSummary(maf_never)$total
smoker_tmb <- getSampleSummary(maf_smoker)$total
never_gene <- getGeneSummary(maf_never)
smoker_gene <- getGeneSummary(maf_smoker)

get_gene_pct <- function(gene, gs, denom) {
  val <- gs %>% filter(Hugo_Symbol == gene) %>% pull(MutatedSamples)
  safe_pct(val, denom)
}

summary_table <- data.frame(
  Feature = c("Samples (n)", "Median mutations per sample", "EGFR mutation frequency",
              "KRAS mutation frequency", "COSMIC_4 mean contribution"),
  Never_Smoker = c(length(never_maf_bc), round(median(never_tmb), 1),
                   paste0(round(get_gene_pct("EGFR", never_gene, length(never_maf_bc)), 1), "%"),
                   paste0(round(get_gene_pct("KRAS", never_gene, length(never_maf_bc)), 1), "%"),
                   ifelse(is.na(c4_mean_never), "NA", paste0(round(c4_mean_never, 1), "%"))),
  Smoker = c(length(smoker_maf_bc), round(median(smoker_tmb), 1),
             paste0(round(get_gene_pct("EGFR", smoker_gene, length(smoker_maf_bc)), 1), "%"),
             paste0(round(get_gene_pct("KRAS", smoker_gene, length(smoker_maf_bc)), 1), "%"),
             ifelse(is.na(c4_mean_smoker), "NA", paste0(round(c4_mean_smoker, 1), "%"))),
  P_value = c("—", formatC(wilcox.test(never_tmb, smoker_tmb)$p.value, format = "e", digits = 2),
              "see mafCompare", "see mafCompare", "see Fig4d")
)

write.csv(summary_table, file.path(DATA, "summary_statistics.csv"), row.names = FALSE, quote = FALSE)
png(file.path(FIGURES, "Summary_Statistics_Table.png"), width = 1600, height = 600, res = 200)
gridExtra::grid.table(summary_table, rows = NULL,
                      theme = gridExtra::ttheme_minimal(
                        core = list(fg_params = list(cex = 0.9)),
                        colhead = list(fg_params = list(cex = 1.0, fontface = "bold"))))
dev.off()
message("Summary table saved.")

# ============================================================
# Final summary
# ============================================================
message("============================================")
message("MAF analysis complete. Figures generated:")
print(list.files(FIGURES, pattern = "^(MAF_|Summary_).*\\.png$"))
message("============================================")
message("Key statistics for presentation:")
message("Never-smoker median mutations: ", median(never_tmb))
message("Smoker median mutations: ", median(smoker_tmb))
message("Top gene in never-smokers: ", never_gene$Hugo_Symbol[1], " (",
        round(never_gene$MutatedSamples[1] / length(never_maf_bc) * 100, 1), "%)")
