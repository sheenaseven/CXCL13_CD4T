suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(harmony)
  library(ComplexHeatmap)
})

## -------------------------------
## 1. Preprocessing and clustering
## -------------------------------

DefaultAssay(merged_seurat_CD4T) <- "RNA"

merged_seurat_CD4T <- merged_seurat_CD4T %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.), npcs = 50) %>%
  RunHarmony(group.by.vars = "orig.ident") %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

merged_seurat_CD4T$cluster <- paste0("sub", merged_seurat_CD4T$RNA_snn_res.0.5)

cluster_levels <- paste0("sub", sort(as.numeric(gsub("sub", "", unique(merged_seurat_CD4T$cluster)))))
merged_seurat_CD4T$cluster <- factor(merged_seurat_CD4T$cluster, levels = cluster_levels)

p_umap <- DimPlot(
  merged_seurat_CD4T,
  reduction = "umap",
  group.by = "cluster",
  label = TRUE,
  raster = FALSE,
  pt.size = 0.1
)

p_umap

## -------------------------------
## 2. Marker dot plot
## -------------------------------

marker_genes <- unique(c(
  "CD3D", "CD4", "CCR7", "LEF1", "IL7R", "GPR183",
  "TNFAIP3", "HSPA1A",
  "FOXP3", "CTLA4",
  "CXCL13", "TOX2", "PDCD1", "CXCR5",
  "GNG4", "TNFRSF18", "ZNRF1", "TNFRSF4",
  "ICOS", "IL2RA", "IL21", "BCL6", "SH2D1A",
  "CD200", "IL13", "IL17RB"
))

marker_genes <- marker_genes[marker_genes %in% rownames(merged_seurat_CD4T)]

p_markers <- DotPlot(
  merged_seurat_CD4T,
  features = marker_genes,
  assay = "RNA",
  group.by = "cluster"
) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p_markers

## -------------------------------
## 3. Cluster proportion by sample
## -------------------------------

df_prop <- as.data.frame(table(
  cluster = merged_seurat_CD4T$cluster,
  sampleID = merged_seurat_CD4T$sampleID
)) %>%
  group_by(sampleID) %>%
  mutate(prop = Freq / sum(Freq) * 100) %>%
  ungroup()

p_prop <- ggplot(df_prop, aes(x = sampleID, y = prop, fill = cluster)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = color_c) +
  labs(x = "CCM sample", y = "Proportion (%)", fill = "Cluster") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_prop

## -------------------------------
## 4. Enrichment analysis
## -------------------------------

contingency_table <- table(
  merged_seurat_CD4T$cluster,
  merged_seurat_CD4T$group
)

chisq_res <- chisq.test(contingency_table)
ro_e <- chisq_res$observed / chisq_res$expected

convert_to_symbol <- function(x) {
  dplyr::case_when(
    x > 1 ~ "++",
    x > 0.8 & x <= 1 ~ "+",
    x >= 0.2 & x <= 0.8 ~ "+",
    x > 0 & x < 0.2 ~ "+/-",
    TRUE ~ "-"
  )
}

mat_symbols <- apply(ro_e, c(1, 2), convert_to_symbol)

p_enrich <- ComplexHeatmap::pheatmap(
  ro_e,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  name = "Ro/e",
  border_color = NA,
  display_numbers = mat_symbols,
  number_color = "black",
  fontsize_number = 10
)

p_enrich

## -------------------------------
## 5. Split sub8 by CXCL13 expression
## -------------------------------

sub8 <- subset(merged_seurat_CD4T, subset = cluster == "sub8")

cxcl13_expr <- GetAssayData(sub8, assay = "RNA", layer = "data")["CXCL13", ]

sub8$newcluster <- ifelse(
  cxcl13_expr > 0,
  "sub8_CXCL13+",
  "sub8_CXCL13-"
)

sub8$newcluster <- factor(
  sub8$newcluster,
  levels = c("sub8_CXCL13-", "sub8_CXCL13+")
)

genes_sub8 <- c(
  "CXCL13", "TOX2", "PDCD1", "CXCR5",
  "TNFRSF18", "ZNRF1", "TNFRSF4", "IL21"
)

genes_sub8 <- genes_sub8[genes_sub8 %in% rownames(sub8)]

p_sub8_dot <- DotPlot(
  sub8,
  features = genes_sub8,
  assay = "RNA",
  group.by = "newcluster"
) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p_sub8_dot

## -------------------------------
## 6. Percent expression and average expression plot
## -------------------------------

plot_df <- p_sub8_dot$data %>%
  filter(features.plot %in% genes_sub8) %>%
  mutate(
    group = factor(id, levels = c("sub8_CXCL13-", "sub8_CXCL13+")),
    pct.exp.raw = pct.exp,
    pct.exp = pct.exp / 100
  )

p_sub8_expr_summary <- ggplot(plot_df, aes(x = group)) +
  geom_col(aes(y = pct.exp), fill = "#F03B20", width = 0.6) +
  geom_col(aes(y = -avg.exp), fill = "#2C7FB8", width = 0.6) +
  geom_text(
    aes(y = pct.exp, label = sprintf("%.1f%%", pct.exp.raw)),
    vjust = 1.3,
    size = 3
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  facet_wrap(~features.plot, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold")
  )

p_sub8_expr_summary

## -------------------------------
## 7. Differential expression
## -------------------------------

Idents(sub8) <- "group"

markers_sub8 <- FindMarkers(
  sub8,
  ident.1 = "CCM",
  test.use = "wilcox",
  logfc.threshold = 0.1,
  min.pct = 0.05
)

markers_sub8 <- markers_sub8 %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    log10_padj = -log10(p_val_adj),
    log10_padj = ifelse(
      is.infinite(log10_padj),
      max(log10_padj[is.finite(log10_padj)], na.rm = TRUE) + 10,
      log10_padj
    )
  )

p_volcano <- ggplot(markers_sub8, aes(x = avg_log2FC, y = log10_padj)) +
  geom_point(aes(color = direction), size = 1, alpha = 0.8) +
  labs(
    x = "Average log2 fold change",
    y = expression(-log[10]("adjusted P value")),
    color = NULL
  ) +
  theme_classic(base_size = 16)

p_volcano