# ============================================================
# TCGA LUAD — Never-Smoker vs Smoker Analysis
# Sherlock-Lung Presentation
# Abid Al Reza, PhD
# ============================================================

# ---- PATHS ----
BASE    <- "/media/wrath/bioinfor_learning/sherlock_lung"
DATA    <- file.path(BASE, "data/TCGA_LUAD")
FIGURES <- file.path(BASE, "figures")
dir.create(FIGURES, showWarnings=FALSE, recursive=TRUE)

# ---- LIBRARIES ----
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggrepel)

# ---- 1. LOAD DATA ----
message("Loading expression data...")
expr_raw <- fread(
  file.path(DATA, "TCGA_LUAD_HiSeqV2.gz"),
  data.table=FALSE
)
rownames(expr_raw) <- expr_raw[,1]
expr_raw <- expr_raw[,-1]
message(paste("Expression:", nrow(expr_raw), "genes x", ncol(expr_raw), "samples"))

message("Loading clinical data...")
clin <- fread(
  file.path(DATA, "TCGA_LUAD_clinical.tsv"),
  data.table=FALSE
)
rownames(clin) <- clin[,1]
message(paste("Clinical:", nrow(clin), "patients"))

# ---- 2. CHECK SMOKING VARIABLE ----
message("Smoking status distribution:")
print(table(clin$tobacco_smoking_history, useNA="always"))

# TCGA tobacco_smoking_history codes:
# 1 = Lifelong non-smoker (<100 cigarettes in lifetime)
# 2 = Current smoker
# 3 = Current reformed smoker >15 years
# 4 = Current reformed smoker ≤15 years
# 5 = Current reformed smoker, duration unknown

# ---- 3. DEFINE GROUPS ----
never_ids <- substr(rownames(clin), 1, 12)[clin$tobacco_smoking_history == 1 & 
                             !is.na(clin$tobacco_smoking_history)]
smoker_ids <- substr(rownames(clin), 1, 12)[clin$tobacco_smoking_history %in% c(2,3,4) & 
                              !is.na(clin$tobacco_smoking_history)]

message(paste("Never-smokers:", length(never_ids)))
message(paste("Smokers:", length(smoker_ids)))

# Match sample IDs between expression and clinical
# Expression samples are TCGA barcodes — take first 12 chars to match patient ID
expr_patient_ids <- substr(colnames(expr_raw), 1, 12)

never_cols <- which(expr_patient_ids %in% never_ids)
smoker_cols <- which(expr_patient_ids %in% smoker_ids)

message(paste("Never-smoker samples in expression:", length(never_cols)))
message(paste("Smoker samples in expression:", length(smoker_cols)))

# Subset expression
expr_never  <- expr_raw[, never_cols]
expr_smoker <- expr_raw[, smoker_cols]
expr_subset <- cbind(expr_never, expr_smoker)

# ---- 4. DIFFERENTIAL EXPRESSION — DESeq2 ----
message("Running DESeq2...")
# Note: Xena HiSeqV2 is log2(x+1) normalized RPKM
# Convert back to integer-like counts for DESeq2
# OR use limma/voom which handles normalized data

library(limma)

# Create group factor
group <- factor(c(
  rep("Never_Smoker", length(never_cols)),
  rep("Smoker", length(smoker_cols))
))

# Design matrix
design <- model.matrix(~group)

# Fit linear model
fit <- lmFit(expr_subset, design)
fit <- eBayes(fit)

# Get results
de_results <- topTable(fit, 
                        coef=2,
                        number=Inf,
                        sort.by="P")

message(paste("DE genes (FDR<0.05):", sum(de_results$adj.P.Val < 0.05)))
message(paste("Upregulated in never-smokers:", 
              sum(de_results$adj.P.Val < 0.05 & de_results$logFC > 1)))
message(paste("Downregulated in never-smokers:", 
              sum(de_results$adj.P.Val < 0.05 & de_results$logFC < -1)))

# Save DE results
write.csv(de_results, 
          file.path(DATA, "DE_never_vs_smoker.csv"),
          quote=FALSE)

# ---- FIGURE 1 — VOLCANO PLOT ----
message("Generating Figure 1: Volcano plot...")

de_plot <- de_results %>%
  mutate(
    gene = rownames(.),
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up in Never-Smoker",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down in Never-Smoker",
      TRUE ~ "NS"
    ),
    neg_log10_p = -log10(adj.P.Val)
  )

# Top genes to label
top_genes <- de_plot %>%
  filter(significance != "NS") %>%
  arrange(adj.P.Val) %>%
  head(20)

p_volcano <- ggplot(de_plot, aes(x=logFC, y=neg_log10_p, color=significance)) +
  geom_point(alpha=0.4, size=1) +
  geom_point(data=filter(de_plot, significance != "NS"),
             alpha=0.7, size=1.5) +
  geom_text_repel(data=top_genes,
                  aes(label=gene),
                  size=3, max.overlaps=15,
                  fontface="italic") +
  scale_color_manual(values=c(
    "Up in Never-Smoker"   = "#D6604D",
    "Down in Never-Smoker" = "#4393C3",
    "NS"                   = "grey70"
  )) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  labs(
    title="Differential Gene Expression\nNever-Smoker vs Smoker Lung Adenocarcinoma (TCGA)",
    subtitle=paste0("Never-smokers n=", length(never_cols), 
                    " | Smokers n=", length(smoker_cols),
                    " | FDR<0.05, |logFC|>1"),
    x="log2 Fold Change (Never-Smoker vs Smoker)",
    y="-log10(Adjusted P-value)",
    color=""
  ) +
  theme_classic(base_size=12) +
  theme(
    legend.position="top",
    plot.title=element_text(face="bold"),
    plot.subtitle=element_text(color="grey40")
  )

ggsave(file.path(FIGURES, "Fig1_Volcano_NeverSmoker_vs_Smoker.png"),
       p_volcano, width=10, height=8, dpi=300)
message("Figure 1 saved.")

# ---- 5. IMMUNE DECONVOLUTION — ESTIMATE ----
message("Running ESTIMATE immune deconvolution...")

# Install ESTIMATE if needed
if (!requireNamespace("estimate", quietly=TRUE)) {
  library(utils)
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
library(estimate)

# ESTIMATE requires specific format
# Write expression to temp file
tmp_expr <- file.path(DATA, "expr_for_estimate.gct")

# Convert to GCT format
expr_gct <- cbind(
  NAME=rownames(expr_subset),
  Description=rownames(expr_subset),
  expr_subset
)
write.table(
  rbind(c("#1.2", rep("", ncol(expr_gct)-1)),
        c(nrow(expr_subset), ncol(expr_subset), rep("", ncol(expr_gct)-2)),
        expr_gct),
  tmp_expr, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
)

# Run ESTIMATE
filterCommonGenes(input.f=tmp_expr,
                  output.f=file.path(DATA, "expr_filtered.gct"),
                  id="GeneSymbol")

estimateScore(input.ds=file.path(DATA, "expr_filtered.gct"),
              output.ds=file.path(DATA, "estimate_scores.gct"))

# Load ESTIMATE results
est_scores <- read.table(
  file.path(DATA, "estimate_scores.gct"),
  skip=2,
  header=TRUE,
  sep="\t",
  row.names=1,
  check.names=FALSE
)
est_scores <- est_scores[,-1]
est_scores_t <- as.data.frame(t(est_scores))

# ESTIMATE/read.table may convert TCGA barcodes from hyphen to dot format.
# Normalize back to TCGA-XX-XXXX-XX style before matching clinical IDs.
est_scores_t$sample_id <- gsub("\\.", "-", rownames(est_scores_t))
est_scores_t$sample_id <- sub("^X", "", est_scores_t$sample_id)

# Add smoking status
est_scores_t$patient_id <- substr(est_scores_t$sample_id, 1, 12)

est_scores_t$smoking <- dplyr::case_when(
  est_scores_t$patient_id %in% never_ids  ~ "Never-Smoker",
  est_scores_t$patient_id %in% smoker_ids ~ "Smoker",
  TRUE ~ NA_character_
)

est_scores_t <- est_scores_t %>%
  dplyr::filter(!is.na(smoking))
  
message("ESTIMATE smoking groups after matching:")
print(table(est_scores_t$smoking, useNA="always"))

if (nrow(est_scores_t) == 0) {
  stop("No ESTIMATE samples matched smoking groups. Check sample_id format.")
}

library(ggpubr)
# ---- FIGURE 2 — IMMUNE/STROMAL SCORES ----
message("Generating Figure 2: ESTIMATE immune scores...")

scores_long <- est_scores_t %>%
  select(sample_id, smoking, ImmuneScore, StromalScore, ESTIMATEScore) %>%
  pivot_longer(cols=c(ImmuneScore, StromalScore, ESTIMATEScore),
               names_to="Score_Type",
               values_to="Score")

p_estimate <- ggplot(scores_long, 
                     aes(x=smoking, y=Score, fill=smoking)) +
  geom_violin(alpha=0.7, trim=FALSE) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
  facet_wrap(~Score_Type, scales="free_y") +
  scale_fill_manual(values=c("Never-Smoker"="#D6604D", "Smoker"="#4393C3")) +
  stat_compare_means(method="wilcox.test",
                     label="p.signif",
                     comparisons=list(c("Never-Smoker", "Smoker"))) +
  labs(
    title="Tumor Microenvironment Immune Landscape\nNever-Smoker vs Smoker LUAD (TCGA, ESTIMATE)",
    x="",
    y="ESTIMATE Score",
    fill=""
  ) +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        plot.title=element_text(face="bold"),
        strip.background=element_rect(fill="grey90"))

ggsave(file.path(FIGURES, "Fig2_ESTIMATE_Immune_Scores.png"),
       p_estimate, width=10, height=6, dpi=300)
message("Figure 2 saved.")

# ---- 6. GSEA PATHWAY ANALYSIS ----
message("Running GSEA...")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# Ranked gene list for GSEA
gene_list <- de_results$logFC
names(gene_list) <- rownames(de_results)
gene_list <- sort(gene_list, decreasing=TRUE)

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(names(gene_list),
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Hs.eg.db)

gene_list_entrez <- gene_list[gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID

# GSEA with HALLMARK gene sets
gsea_hallmark <- gseKEGG(
  geneList=gene_list_entrez,
  organism="hsa",
  minGSSize=15,
  maxGSSize=500,
  pvalueCutoff=0.05,
  pAdjustMethod="BH",
  verbose=FALSE
)

# ---- FIGURE 3 — GSEA DOTPLOT ----
message("Generating Figure 3: GSEA...")

if (nrow(gsea_hallmark@result) > 0) {
  p_gsea <- dotplot(gsea_hallmark, 
                    showCategory=20,
                    split=".sign") +
    facet_grid(.~.sign) +
    labs(title="KEGG Pathway Enrichment\nNever-Smoker vs Smoker LUAD",
         subtitle="Gene Set Enrichment Analysis (GSEA)") +
    theme(plot.title=element_text(face="bold"))
  
  ggsave(file.path(FIGURES, "Fig3_GSEA_Pathways.png"),
         p_gsea, width=14, height=10, dpi=300)
  message("Figure 3 saved.")
} else {
  message("No significant GSEA pathways found — try MSigDB Hallmark sets")
}

# ---- 7. SURVIVAL ANALYSIS ----
message("Running survival analysis...")
library(survival)
library(survminer)

# Merge immune score with survival data
surv_data <- est_scores_t %>%
  mutate(patient_id=substr(sample_id, 1, 12)) %>%
  left_join(
    clin %>% 
      dplyr::select(
        sampleID, 
        OS.time=days_to_death,
        OS.time2=days_to_last_followup,
        OS.status=vital_status,
        days_to_last_followup,
        tobacco_smoking_history
      ) %>%
      mutate(patient_id=substr(sampleID, 1, 12)),
    by="patient_id"
  ) %>%
  filter((!is.na(OS.time) | !is.na(OS.time2)), !is.na(ImmuneScore)) %>%
  mutate(
    OS.time = as.numeric(ifelse(!is.na(OS.time), OS.time, OS.time2)) / 365,
    OS.status = ifelse(!is.na(OS.status) & OS.status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS.time), OS.time > 0) %>%
  mutate(
    ImmuneHigh = ifelse(
      ImmuneScore > median(ImmuneScore, na.rm=TRUE),
      "Immune High",
      "Immune Low"
    )
  )

# Kaplan-Meier by immune score
km_fit <- survfit(Surv(OS.time, OS.status) ~ ImmuneHigh,
                  data=surv_data)

# ---- FIGURE 4 — SURVIVAL BY IMMUNE SCORE ----
message("Generating Figure 4: Survival analysis...")

p_surv <- ggsurvplot(
  km_fit,
  data=surv_data,
  pval=TRUE,
  conf.int=TRUE,
  risk.table=TRUE,
  palette=c("#D6604D", "#4393C3"),
  title="Overall Survival by Tumor Immune Infiltration\nLUAD (TCGA)",
  xlab="Time (Years)",
  ylab="Overall Survival Probability",
  legend.title="",
  legend.labs=c("Immune High", "Immune Low"),
  ggtheme=theme_classic(base_size=12)
)

png(file.path(FIGURES, "Fig4_Survival_ImmuneScore.png"), width=3000, height=2400, res=300)
print(p_surv)
dev.off()
message("Figure 4 saved.")

# ---- FINAL SUMMARY ----
message("==============================================")
message("Analysis complete. Figures saved:")
list.files(FIGURES)
message("==============================================")
