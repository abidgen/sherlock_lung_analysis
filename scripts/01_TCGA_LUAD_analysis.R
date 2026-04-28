# ============================================================
# TCGA LUAD — Never-Smoker vs Smoker Expression Analysis
# Sherlock-Lung Presentation
# Abid Al Reza, PhD
# ============================================================

# ---- Paths ----
BASE <- Sys.getenv("SHERLOCK_LUNG_BASE", unset = getwd())
DATA <- file.path(BASE, "data", "TCGA_LUAD")
FIGURES <- file.path(BASE, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----
required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggrepel", "ggpubr",
  "limma", "estimate", "clusterProfiler", "enrichplot", "org.Hs.eg.db",
  "survival", "survminer"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them in the conda environment and rerun.")
}

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# ---- Helper functions ----
check_file <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) stop("Missing or empty file: ", path)
  invisible(path)
}

patient12 <- function(x) substr(gsub("\\.", "-", x), 1, 12)

# ---- Input files ----
expr_file <- check_file(file.path(DATA, "TCGA_LUAD_HiSeqV2.gz"))
clin_file <- check_file(file.path(DATA, "TCGA_LUAD_clinical.tsv"))

# ============================================================
# 1. Load expression and clinical data
# ============================================================
message("Loading expression data...")
expr_raw <- fread(expr_file, data.table = FALSE)
rownames(expr_raw) <- expr_raw[[1]]
expr_raw <- expr_raw[, -1, drop = FALSE]
message("Expression: ", nrow(expr_raw), " genes x ", ncol(expr_raw), " samples")

message("Loading clinical data...")
clin <- fread(clin_file, data.table = FALSE)
if (!"sampleID" %in% colnames(clin)) colnames(clin)[1] <- "sampleID"
clin <- clin %>% mutate(patient_id = patient12(sampleID))
message("Clinical: ", nrow(clin), " rows")

# ============================================================
# 2. Define smoking groups and match expression samples
# ============================================================
message("Smoking status distribution:")
print(table(clin$tobacco_smoking_history, useNA = "always"))

clin_smoking <- clin %>%
  mutate(smoking = case_when(
    tobacco_smoking_history == 1 ~ "Never-Smoker",
    tobacco_smoking_history %in% c(2, 3, 4) ~ "Smoker",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(smoking)) %>%
  distinct(patient_id, .keep_all = TRUE)

never_ids <- clin_smoking %>% filter(smoking == "Never-Smoker") %>% pull(patient_id)
smoker_ids <- clin_smoking %>% filter(smoking == "Smoker") %>% pull(patient_id)

expr_patient_ids <- patient12(colnames(expr_raw))
never_cols <- which(expr_patient_ids %in% never_ids)
smoker_cols <- which(expr_patient_ids %in% smoker_ids)

message("Never-smoker samples in expression: ", length(never_cols))
message("Smoker samples in expression: ", length(smoker_cols))
if (length(never_cols) < 2 || length(smoker_cols) < 2) stop("Too few matched expression samples.")

expr_subset <- cbind(expr_raw[, never_cols, drop = FALSE],
                     expr_raw[, smoker_cols, drop = FALSE])
group <- factor(c(rep("Never_Smoker", length(never_cols)),
                  rep("Smoker", length(smoker_cols))))

# ============================================================
# 3. Differential expression with limma
# ============================================================
message("Running limma differential expression...")
library(limma)

design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr_subset, design))
de_results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

message("DE genes FDR<0.05: ", sum(de_results$adj.P.Val < 0.05, na.rm = TRUE))
message("Up in never-smokers: ", sum(de_results$adj.P.Val < 0.05 & de_results$logFC > 1, na.rm = TRUE))
message("Down in never-smokers: ", sum(de_results$adj.P.Val < 0.05 & de_results$logFC < -1, na.rm = TRUE))
write.csv(de_results, file.path(DATA, "DE_never_vs_smoker.csv"), quote = FALSE)

# ---- Figure 1: volcano ----
message("Generating Figure 1: Volcano plot...")
de_plot <- de_results %>%
  mutate(
    gene = rownames(.),
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up in Never-Smoker",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down in Never-Smoker",
      TRUE ~ "NS"
    ),
    neg_log10_p = -log10(adj.P.Val)
  )

top_genes <- de_plot %>% filter(significance != "NS") %>% arrange(adj.P.Val) %>% head(20)

p_volcano <- ggplot(de_plot, aes(x = logFC, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_point(data = filter(de_plot, significance != "NS"), alpha = 0.7, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3,
                  max.overlaps = 15, fontface = "italic") +
  scale_color_manual(values = c("Up in Never-Smoker" = "#D6604D",
                                "Down in Never-Smoker" = "#4393C3",
                                "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(title = "Differential Gene Expression\nNever-Smoker vs Smoker LUAD (TCGA)",
       subtitle = paste0("Never-smokers n=", length(never_cols),
                         " | Smokers n=", length(smoker_cols),
                         " | FDR<0.05, |logFC|>1"),
       x = "log2 Fold Change (Never-Smoker vs Smoker)",
       y = "-log10(Adjusted P-value)", color = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))

ggsave(file.path(FIGURES, "Fig1_Volcano_NeverSmoker_vs_Smoker.png"),
       p_volcano, width = 10, height = 8, dpi = 200)
message("Figure 1 saved.")

# ============================================================
# 4. ESTIMATE immune/stromal scores
# ============================================================
message("Running ESTIMATE immune deconvolution...")
library(estimate)

tmp_expr <- file.path(DATA, "expr_for_estimate.gct")

expr_gct <- cbind(
  NAME = rownames(expr_subset),
  Description = rownames(expr_subset),
  expr_subset
)

write.table(
  rbind(
    c("#1.2", rep("", ncol(expr_gct) - 1)),
    c(nrow(expr_subset), ncol(expr_subset), rep("", ncol(expr_gct) - 2)),
    expr_gct
  ),
  tmp_expr,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

filterCommonGenes(
  input.f = tmp_expr,
  output.f = file.path(DATA, "expr_filtered.gct"),
  id = "GeneSymbol"
)

estimateScore(
  input.ds = file.path(DATA, "expr_filtered.gct"),
  output.ds = file.path(DATA, "estimate_scores.gct")
)

est_scores <- read.table(
  file.path(DATA, "estimate_scores.gct"),
  skip = 2,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

est_scores <- est_scores[, -1, drop = FALSE]
est_scores_t <- as.data.frame(t(est_scores))

est_scores_t$sample_id <- gsub("\\.", "-", rownames(est_scores_t))
est_scores_t$sample_id <- sub("^X", "", est_scores_t$sample_id)
est_scores_t$patient_id <- patient12(est_scores_t$sample_id)

est_scores_t <- est_scores_t %>%
  mutate(smoking = case_when(
    patient_id %in% never_ids ~ "Never-Smoker",
    patient_id %in% smoker_ids ~ "Smoker",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(smoking))

message("ESTIMATE matched samples: ", nrow(est_scores_t))
print(table(est_scores_t$smoking, useNA = "always"))

if (nrow(est_scores_t) == 0) {
  stop("No ESTIMATE samples matched smoking groups. Check sample_id format.")
}

# ---- Figure 2: ESTIMATE scores ----
message("Generating Figure 2: ESTIMATE immune scores...")

scores_long <- est_scores_t %>%
  dplyr::select(sample_id, smoking, ImmuneScore, StromalScore, ESTIMATEScore) %>%
  tidyr::pivot_longer(
    cols = c(ImmuneScore, StromalScore, ESTIMATEScore),
    names_to = "Score_Type",
    values_to = "Score"
  )

p_estimate <- ggplot(scores_long, aes(x = smoking, y = Score, fill = smoking)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ Score_Type, scales = "free_y") +
  scale_fill_manual(values = c("Never-Smoker" = "#D6604D", "Smoker" = "#4393C3")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("Never-Smoker", "Smoker"))
  ) +
  labs(
    title = "Tumor Microenvironment Immune Landscape\nNever-Smoker vs Smoker LUAD (TCGA, ESTIMATE)",
    x = "",
    y = "ESTIMATE Score",
    fill = ""
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  file.path(FIGURES, "Fig2_ESTIMATE_Immune_Scores.png"),
  p_estimate,
  width = 10,
  height = 6,
  dpi = 200
)

message("Figure 2 saved.")

# ============================================================
# 5. GSEA pathway analysis
# ============================================================
message("Running GSEA...")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

gene_list <- de_results$logFC
names(gene_list) <- rownames(de_results)
gene_list <- sort(gene_list, decreasing = TRUE)

gene_ids <- bitr(names(gene_list), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)

gene_df <- data.frame(SYMBOL = names(gene_list), logFC = as.numeric(gene_list)) %>%
  inner_join(gene_ids, by = "SYMBOL") %>%
  filter(!is.na(ENTREZID), !is.na(logFC)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_entrez <- gene_df$logFC
names(gene_list_entrez) <- as.character(gene_df$ENTREZID)
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

if (requireNamespace("msigdbr", quietly = TRUE)) {
  message("Using MSigDB Hallmark gene sets via msigdbr.")
  library(msigdbr)

  msig_h <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )

  hallmark_sets <- msig_h %>%
    dplyr::select(gs_name, gene = ncbi_gene) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::mutate(gene = as.character(gene))

  gsea_res <- GSEA(
    geneList = gene_list_entrez,
    TERM2GENE = hallmark_sets,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  fig3_name <- "Fig3_GSEA_Hallmark.png"
  fig3_title <- "Hallmark Pathway Enrichment\nNever-Smoker vs Smoker LUAD (TCGA)"
} else {
  message("Package msigdbr not found; using KEGG GSEA fallback.")
  gsea_res <- gseKEGG(geneList = gene_list_entrez, organism = "hsa",
                      minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05,
                      pAdjustMethod = "BH", verbose = FALSE)
  fig3_name <- "Fig3_GSEA_KEGG.png"
  fig3_title <- "KEGG Pathway Enrichment\nNever-Smoker vs Smoker LUAD (TCGA)"
}

message("Significant pathways: ", nrow(gsea_res@result))

# ---- Figure 3: GSEA dotplot ----
message("Generating Figure 3: GSEA...")
if (nrow(gsea_res@result) > 0) {
  p_gsea <- dotplot(gsea_res, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = fig3_title, subtitle = "Gene Set Enrichment Analysis") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 9),
          strip.background = element_rect(fill = "grey90"))
  ggsave(file.path(FIGURES, fig3_name), p_gsea, width = 14, height = 10, dpi = 200)
  message("Figure 3 saved.")
} else {
  message("No significant pathways at FDR<0.05; skipping GSEA figure.")
}

# ============================================================
# 6. Survival by immune score
# ============================================================
message("Running survival analysis...")
library(survival)
library(survminer)

clin_surv <- clin %>%
  dplyr::select(
    patient_id,
    sampleID,
    OS.time = days_to_death,
    OS.time2 = days_to_last_followup,
    OS.status = vital_status,
    tobacco_smoking_history
  ) %>%
  dplyr::distinct(patient_id, .keep_all = TRUE)

surv_data <- est_scores_t %>%
  dplyr::left_join(clin_surv, by = "patient_id") %>%
  dplyr::filter((!is.na(OS.time) | !is.na(OS.time2)), !is.na(ImmuneScore)) %>%
  dplyr::mutate(
    OS.time = suppressWarnings(as.numeric(ifelse(!is.na(OS.time), OS.time, OS.time2))) / 365,
    OS.status = ifelse(!is.na(OS.status) & OS.status == "DECEASED", 1, 0),
    ImmuneHigh = ifelse(
      ImmuneScore > median(ImmuneScore, na.rm = TRUE),
      "Immune High",
      "Immune Low"
    ),
    smoking_grp = dplyr::case_when(
      tobacco_smoking_history == 1 ~ "Never_Smoker",
      tobacco_smoking_history %in% c(2, 3, 4) ~ "Smoker",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(OS.time), OS.time > 0)

message("Survival analysis samples: ", nrow(surv_data))
message("Events (deaths): ", sum(surv_data$OS.status, na.rm = TRUE))

if (nrow(surv_data) > 0 && sum(surv_data$OS.status, na.rm = TRUE) >= 5) {
  km_fit <- survfit(Surv(OS.time, OS.status) ~ ImmuneHigh, data = surv_data)
  p_surv <- ggsurvplot(km_fit, data = surv_data,
                       pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                       palette = c("#D6604D", "#4393C3"),
                       title = "Overall Survival by Tumor Immune Infiltration\nLUAD (TCGA)",
                       xlab = "Time (Years)", ylab = "Overall Survival Probability",
                       legend.title = "", legend.labs = c("Immune High", "Immune Low"),
                       ggtheme = theme_classic(base_size = 12))
  png(file.path(FIGURES, "Fig4_Survival_ImmuneScore.png"), width = 3000, height = 2400, res = 300)
  print(p_surv)
  dev.off()
  message("Figure 4 saved.")
} else {
  message("Too few survival observations/events; skipping Fig4.")
}

message("==============================================")
message("Expression analysis complete. Figures saved:")
print(list.files(FIGURES, pattern = "^Fig[0-9].*\\.png$"))
message("==============================================")
