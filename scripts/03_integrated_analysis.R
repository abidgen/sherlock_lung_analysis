# ============================================================
# TCGA LUAD — Integrated Driver-Expression and CNA Analysis
# EGFR-mutant vs wildtype | CNA support
# Abid Al Reza, PhD — Sherlock-Lung Presentation
# ============================================================

# ---- Paths ----
BASE <- Sys.getenv("SHERLOCK_LUNG_BASE", unset = getwd())
DATA <- file.path(BASE, "data", "TCGA_LUAD")
FIGURES <- file.path(BASE, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----
required_pkgs <- c("maftools", "data.table", "dplyr", "ggplot2", "ggrepel", "ggpubr", "patchwork", "limma", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them in the conda environment and rerun.")
}

library(maftools)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(limma)
library(tidyr)

# ---- Helpers ----
check_file <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) stop("Missing or empty file: ", path)
  invisible(path)
}
patient12 <- function(x) substr(gsub("\\.", "-", x), 1, 12)
read_cna <- function(data_dir) {
  candidates <- c(
    file.path(data_dir, "TCGA_LUAD_CNA.tsv"),
    file.path(data_dir, "TCGA_LUAD_CNA.tsv.gz"),
    file.path(data_dir, "TCGA_PANCAN_CNA.tsv.gz")
  )
  hit <- candidates[file.exists(candidates) & file.size(candidates) > 100000]
  if (length(hit) == 0) return(NULL)
  path <- hit[1]
  message("Loading CNA file: ", basename(path))
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("zcat", shQuote(path)), data.table = FALSE)
  } else {
    fread(path, data.table = FALSE)
  }
}

# ---- Input files ----
clin_file <- check_file(file.path(DATA, "TCGA_LUAD_clinical.tsv"))
expr_file <- check_file(file.path(DATA, "TCGA_LUAD_HiSeqV2.gz"))
maf_file <- check_file(file.path(DATA, "TCGA_LUAD_somatic.maf.gz"))

# ============================================================
# 1. Load shared data
# ============================================================
message("Loading shared data...")
clin <- read.delim(clin_file, header = TRUE, sep = "\t")
if (!"sampleID" %in% colnames(clin)) colnames(clin)[1] <- "sampleID"
clin <- clin %>% mutate(patient_id = patient12(sampleID))

expr_raw <- fread(expr_file, data.table = FALSE)
rownames(expr_raw) <- expr_raw[[1]]
expr_raw <- expr_raw[, -1, drop = FALSE]
expr_patient_ids <- patient12(colnames(expr_raw))

luad_maf <- read.maf(maf = maf_file, verbose = FALSE)

never_ids <- clin$patient_id[clin$tobacco_smoking_history == 1 & !is.na(clin$tobacco_smoking_history)]
smoker_ids <- clin$patient_id[clin$tobacco_smoking_history %in% c(2, 3, 4) & !is.na(clin$tobacco_smoking_history)]

get_mutated_patients <- function(maf, gene) {
  maf@data %>%
    filter(Hugo_Symbol == gene) %>%
    mutate(patient_id = patient12(Tumor_Sample_Barcode)) %>%
    pull(patient_id) %>%
    unique()
}

egfr_mutated <- get_mutated_patients(luad_maf, "EGFR")
kras_mutated <- get_mutated_patients(luad_maf, "KRAS")
tp53_mutated <- get_mutated_patients(luad_maf, "TP53")

sample_info <- data.frame(patient_id = unique(expr_patient_ids), stringsAsFactors = FALSE) %>%
  left_join(clin %>% select(patient_id, tobacco_smoking_history, sampleID) %>% distinct(patient_id, .keep_all = TRUE),
            by = "patient_id") %>%
  mutate(
    smoking = case_when(
      tobacco_smoking_history == 1 ~ "Never_Smoker",
      tobacco_smoking_history %in% c(2, 3, 4) ~ "Smoker",
      TRUE ~ NA_character_
    ),
    EGFR_status = ifelse(patient_id %in% egfr_mutated, "EGFR_mutant", "EGFR_wildtype"),
    KRAS_status = ifelse(patient_id %in% kras_mutated, "KRAS_mutant", "KRAS_wildtype"),
    TP53_status = ifelse(patient_id %in% tp53_mutated, "TP53_mutant", "TP53_wildtype")
  ) %>%
  filter(!is.na(smoking))

matched_idx <- match(sample_info$patient_id, expr_patient_ids)
valid <- !is.na(matched_idx)
expr_matched <- expr_raw[, matched_idx[valid], drop = FALSE]
sample_info <- sample_info[valid, ]

message("Integrated expression dataset: ", nrow(sample_info), " samples")
message("EGFR mutant: ", sum(sample_info$EGFR_status == "EGFR_mutant"))
message("EGFR wildtype: ", sum(sample_info$EGFR_status == "EGFR_wildtype"))

# ============================================================
# 2. EGFR expression by mutation status
# ============================================================
message("Generating EGFR expression by mutation status...")
if (!"EGFR" %in% rownames(expr_matched)) stop("EGFR row not found in expression matrix.")

egfr_expr <- data.frame(
  patient_id = sample_info$patient_id,
  EGFR_expr = as.numeric(expr_matched["EGFR", ]),
  EGFR_status = sample_info$EGFR_status,
  KRAS_status = sample_info$KRAS_status,
  smoking = sample_info$smoking,
  stringsAsFactors = FALSE
)

p_egfr_expr <- ggplot(egfr_expr, aes(x = EGFR_status, y = EGFR_expr, fill = EGFR_status)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ smoking) +
  scale_fill_manual(values = c("EGFR_mutant" = "#D6604D", "EGFR_wildtype" = "#4393C3")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(title = "EGFR Expression by Mutation Status",
       subtitle = "Stratified by smoking history — TCGA LUAD",
       x = "", y = "EGFR Expression (log2 RPKM)", fill = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(FIGURES, "Fig5_EGFR_Expression_by_Mutation.png"),
       p_egfr_expr, width = 10, height = 6, dpi = 200)
message("Figure 5 saved.")

# ============================================================
# 3. EGFR-mutant vs wildtype DE in never-smokers
# ============================================================
message("Running EGFR-mutant vs wildtype DE in never-smokers...")
never_mask <- sample_info$smoking == "Never_Smoker"
expr_never_int <- expr_matched[, never_mask, drop = FALSE]
sample_never_int <- sample_info[never_mask, ]
n_mut <- sum(sample_never_int$EGFR_status == "EGFR_mutant")
n_wt <- sum(sample_never_int$EGFR_status == "EGFR_wildtype")
message("Never-smoker EGFR: mutant = ", n_mut, " | wildtype = ", n_wt)

if (n_mut >= 5 && n_wt >= 5) {
  egfr_grp <- factor(sample_never_int$EGFR_status, levels = c("EGFR_wildtype", "EGFR_mutant"))
  fit_egfr <- eBayes(lmFit(expr_never_int, model.matrix(~ egfr_grp)))
  de_egfr <- topTable(fit_egfr, coef = 2, number = Inf, sort.by = "P")
  write.csv(de_egfr, file.path(DATA, "DE_EGFR_mutant_vs_wildtype.csv"), quote = FALSE)
  message("EGFR DE FDR<0.05: ", sum(de_egfr$adj.P.Val < 0.05, na.rm = TRUE))

  de_egfr_plot <- de_egfr %>%
    mutate(
      gene = rownames(.),
      significance = case_when(
        adj.P.Val < 0.05 & logFC > 1 ~ "Up in EGFR-mutant",
        adj.P.Val < 0.05 & logFC < -1 ~ "Down in EGFR-mutant",
        TRUE ~ "NS"
      ),
      neg_log10_p = -log10(adj.P.Val + 1e-10)
    )
  top_egfr <- de_egfr_plot %>% filter(significance != "NS") %>% arrange(adj.P.Val) %>% head(20)

  p_egfr_volcano <- ggplot(de_egfr_plot, aes(x = logFC, y = neg_log10_p, color = significance)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_point(data = filter(de_egfr_plot, significance != "NS"), alpha = 0.8, size = 1.5) +
    geom_text_repel(data = top_egfr, aes(label = gene), size = 3, max.overlaps = 15, fontface = "italic") +
    scale_color_manual(values = c("Up in EGFR-mutant" = "#D6604D",
                                  "Down in EGFR-mutant" = "#4393C3",
                                  "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(title = "EGFR-Mutant vs EGFR-Wildtype Expression\nNever-Smoker LUAD (TCGA)",
         subtitle = paste0("EGFR-mutant n=", n_mut, " | EGFR-wildtype n=", n_wt,
                           " | FDR<0.05, |logFC|>1"),
         x = "log2 Fold Change (EGFR-mutant vs wildtype)",
         y = "-log10(Adjusted P-value)", color = "") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
  ggsave(file.path(FIGURES, "Fig5b_EGFR_Mutant_vs_Wildtype_Volcano.png"),
         p_egfr_volcano, width = 10, height = 8, dpi = 200)
  message("Figure 5b saved.")
} else {
  message("Insufficient EGFR mutant/wildtype samples for never-smoker DE; skipping Fig5b.")
}

# ============================================================
# 4. Copy number alteration analysis
# ============================================================
message("Running CNA analysis...")
cna <- read_cna(DATA)

if (is.null(cna)) {
  message("CNA file not available — skipping CNA analysis.")
  message("Optional file names: TCGA_LUAD_CNA.tsv, TCGA_LUAD_CNA.tsv.gz, or TCGA_PANCAN_CNA.tsv.gz")
} else {
  gene_col <- colnames(cna)[1]
  rownames(cna) <- cna[[gene_col]]
  cna <- cna[, -1, drop = FALSE]
  cna_patient_ids <- patient12(colnames(cna))

  never_cna_cols <- which(cna_patient_ids %in% never_ids)
  smoker_cna_cols <- which(cna_patient_ids %in% smoker_ids)
  message("CNA never-smoker samples: ", length(never_cna_cols))
  message("CNA smoker samples: ", length(smoker_cna_cols))

  if (length(never_cna_cols) >= 5 && length(smoker_cna_cols) >= 5) {
    luad_genes <- c("EGFR", "KRAS", "TP53", "STK11", "KEAP1", "MET", "ERBB2", "RET",
                    "CDKN2A", "RB1", "SMAD4", "NF1", "MYC", "CCND1")
    genes_present <- luad_genes[luad_genes %in% rownames(cna)]

    cna_num <- as.data.frame(lapply(cna[genes_present, c(never_cna_cols, smoker_cna_cols), drop = FALSE],
                                    function(x) suppressWarnings(as.numeric(x))))
    rownames(cna_num) <- genes_present
    n_never <- length(never_cna_cols)
    cna_never <- cna_num[, seq_len(n_never), drop = FALSE]
    cna_smoker <- cna_num[, (n_never + 1):ncol(cna_num), drop = FALSE]

    cna_summary <- data.frame(
      Gene = genes_present,
      Amp_Never = round(rowMeans(cna_never >= 2, na.rm = TRUE) * 100, 1),
      Amp_Smoker = round(rowMeans(cna_smoker >= 2, na.rm = TRUE) * 100, 1),
      Del_Never = round(rowMeans(cna_never <= -2, na.rm = TRUE) * 100, 1),
      Del_Smoker = round(rowMeans(cna_smoker <= -2, na.rm = TRUE) * 100, 1)
    ) %>% arrange(desc(Amp_Never))

    write.csv(cna_summary, file.path(DATA, "CNA_summary.csv"), row.names = FALSE, quote = FALSE)

    cna_amp_long <- cna_summary %>%
      pivot_longer(cols = c(Amp_Never, Amp_Smoker), names_to = "Group", values_to = "Pct") %>%
      mutate(Group = recode(Group, "Amp_Never" = "Never-Smoker", "Amp_Smoker" = "Smoker"))

    p_cna_amp <- ggplot(cna_amp_long, aes(x = reorder(Gene, -Pct), y = Pct, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Never-Smoker" = "#D6604D", "Smoker" = "#4393C3")) +
      labs(title = "Gene Amplification Frequency by Smoking Status\nLUAD (TCGA, GISTIC2)",
           subtitle = "Gene-level thresholded CNA: amplification defined as GISTIC >= 2",
           x = "Gene", y = "Amplification Frequency (%)", fill = "") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            plot.title = element_text(face = "bold"), legend.position = "top")

    cna_del_long <- cna_summary %>%
      pivot_longer(cols = c(Del_Never, Del_Smoker), names_to = "Group", values_to = "Pct") %>%
      mutate(Group = recode(Group, "Del_Never" = "Never-Smoker", "Del_Smoker" = "Smoker"))

    p_cna_del <- ggplot(cna_del_long, aes(x = reorder(Gene, -Pct), y = Pct, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Never-Smoker" = "#D6604D", "Smoker" = "#4393C3")) +
      labs(title = "Deep Deletion Frequency by Smoking Status",
           subtitle = "Gene-level thresholded CNA: deep deletion defined as GISTIC <= -2",
           x = "Gene", y = "Deep Deletion Frequency (%)", fill = "") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            plot.title = element_text(face = "bold"), legend.position = "top")

    combined_cna <- p_cna_amp / p_cna_del +
      plot_annotation(title = "Copy Number Alteration Landscape — TCGA LUAD",
                      theme = theme(plot.title = element_text(face = "bold", size = 14)))
    ggsave(file.path(FIGURES, "Fig7_CNA_Landscape.png"), combined_cna, width = 14, height = 12, dpi = 200)
    message("Figure 7 saved.")

    if ("EGFR" %in% rownames(cna)) {
      egfr_cna_status <- data.frame(
        patient_id = cna_patient_ids,
        EGFR_CNA = suppressWarnings(as.numeric(cna["EGFR", ])),
        stringsAsFactors = FALSE
      ) %>%
        mutate(EGFR_amp = ifelse(EGFR_CNA >= 2, "Amplified", "Not_Amplified"),
               patient_id_12 = patient12(patient_id)) %>%
        left_join(sample_info %>% select(patient_id, smoking, EGFR_status),
                  by = c("patient_id_12" = "patient_id")) %>%
        filter(!is.na(smoking))

      egfr_combined <- egfr_cna_status %>%
        count(smoking, EGFR_status, EGFR_amp) %>%
        group_by(smoking) %>%
        mutate(pct = n / sum(n) * 100)

      p_egfr_cna <- ggplot(egfr_combined,
                            aes(x = smoking, y = pct, fill = interaction(EGFR_status, EGFR_amp))) +
        geom_bar(stat = "identity") +
        scale_fill_manual(
          values = c("EGFR_mutant.Amplified" = "#7B1D1D",
                     "EGFR_mutant.Not_Amplified" = "#D6604D",
                     "EGFR_wildtype.Amplified" = "#08519C",
                     "EGFR_wildtype.Not_Amplified" = "#4393C3"),
          labels = c("Mutant + Amplified", "Mutant only", "Amplified only", "Neither")
        ) +
        labs(title = "EGFR Mutation and Amplification Status by Smoking",
             x = "", y = "% of samples", fill = "EGFR status") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold"), legend.position = "right")
      ggsave(file.path(FIGURES, "Fig7b_EGFR_Mutation_Amplification.png"),
             p_egfr_cna, width = 10, height = 6, dpi = 200)
      message("Figure 7b saved.")
    }
  } else {
    message("Too few CNA samples matched smoking groups; skipping CNA figures.")
  }
}

# ============================================================
# Final summary
# ============================================================
message("==============================================")
message("Integrated analysis complete. Figures:")
print(list.files(FIGURES, pattern = "^Fig[5-9].*\\.png$"))
message("==============================================")
