library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(SCORPIUS)
library(ggplot2)
library(dplyr)

# pbmc_mono <- readRDS("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\sub.RDS")
# print(table(pbmc_mono$celltype))
# meta_sub<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\meta_sub.RDS")
# pbmc_mono@meta.data<-meta_sub
#saveRDS(pbmc_mono,file="E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\pbmc_mono.RDS")
pbmc_mono<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\pbmc_mono.RDS")

######1 slingshot
sce <- as.SingleCellExperiment(pbmc_mono, assay = "RNA")
reducedDim(sce, "UMAP") <- Embeddings(pbmc_mono, "harmony_umap")
sce_sling <- slingshot(sce,
                       reducedDim    = "UMAP",
                       clusterLabels = sce$celltype,
                       approx_points = 1500)

sling_pt <- slingPseudotime(sce_sling)
if (is.matrix(sling_pt)) {
  pbmc_mono$sling_pseudotime <- rowMeans(sling_pt, na.rm = TRUE)
  cat("Lineages:", ncol(sling_pt), "| using mean\n")
} else {
  pbmc_mono$sling_pseudotime <- sling_pt
  cat("Single lineage\n")
}
pbmc_mono$sling_pseudotime[is.nan(pbmc_mono$sling_pseudotime)] <- NA
cat("Slingshot range:", range(pbmc_mono$sling_pseudotime, na.rm = TRUE), "\n")
cat("NA cells:", sum(is.na(pbmc_mono$sling_pseudotime)), "\n")
sum(is.na(pbmc_mono$sling_pseudotime))#0
length(pbmc_mono$sling_pseudotime)#25086
ncol(pbmc_mono)#25086

#rm(sce, sce_sling, sling_pt); gc()


######2 SCORPIUS
library(SCORPIUS)
expr_mat <- as.matrix(t(GetAssayData(pbmc_mono, assay = "RNA", slot = "data")))
space <- reduce_dimensionality(expr_mat, "spearman")
traj  <- infer_trajectory(space)
pbmc_mono$scorpius_pseudotime <- traj$time[colnames(pbmc_mono)]
sum(is.na(pbmc_mono$scorpius_pseudotime))#0
length(pbmc_mono$scorpius_pseudotime)#25086
ncol(pbmc_mono)#25086

######3 Monocle3
library(monocle3)
library(ggplot2)
library(dplyr)
counts <- GetAssayData(pbmc_mono, assay = "RNA", slot = "counts")
meta <- pbmc_mono@meta.data
gene_info <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
cds <- new_cell_data_set(counts, cell_metadata = meta, gene_metadata = gene_info)
cat("CDS:", nrow(cds), "genes x", ncol(cds), "cells\n")

cds <- preprocess_cds(cds, num_dim = 30, verbose = FALSE)
cds <- reduce_dimension(cds, reduction_method = "UMAP", verbose = FALSE)
cds <- cluster_cells(cds, resolution = 1e-5, verbose = FALSE)
cat("Partitions:", length(unique(partitions(cds))), "\n")
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)

cd14_expr <- counts["CD14", ]
root_cell <- names(which.max(cd14_expr))
cat("Root cell (highest CD14):", root_cell, "\n")# max expression cell_id was selected as root

library(igraph)
graph_nodes <- names(V(cds@principal_graph$UMAP))
cat("Graph node examples:", head(graph_nodes), "\n")
cat("Total nodes:", length(graph_nodes), "\n")

root_idx <- as.integer(
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[root_cell, ])
root_node <- graph_nodes[root_idx]
cat("Root vertex index:", root_idx, "| node name:", root_node, "\n")
cds <- order_cells(cds, root_pr_nodes = root_node, verbose = FALSE)

mono_pt <- pseudotime(cds, reduction_method = "UMAP")
pbmc_mono$monocle3_pseudotime <- mono_pt[colnames(pbmc_mono)]
cat("Monocle3 range:", round(range(pbmc_mono$monocle3_pseudotime, na.rm = TRUE), 3), "\n")
cat("NA cells:", sum(is.na(pbmc_mono$monocle3_pseudotime)), "\n")


######4 DiffusionMap
dpt<-read.csv("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\dpt_pseudotime.csv")
ncol(pbmc_mono)
dim(dpt)
pbmc_mono@meta.data$dpt_pseudotime<-dpt[match(rownames(pbmc_mono@meta.data),dpt$X),"dpt_pseudotime"]


######5 plantir
pal<-read.csv("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\palantir_pseudotime.csv")
ncol(pbmc_mono)
dim(pal)
pbmc_mono@meta.data$palantir_pseudotime<-pal[match(rownames(pbmc_mono@meta.data),pal$X),"X0"]
sum(is.na(pbmc_mono$palantir_pseudotime))


######5 check data
check_pseudotime <- function(seurat_obj, pt_col, ct_col = "celltype",
                             colors = c("CD14 Monocyte" = "#6699CC", "CD16 Monocyte" = "#BB0021")) {
  library(ggplot2); library(dplyr)
  pt_df <- FetchData(seurat_obj, vars = c(pt_col, ct_col))
  colnames(pt_df) <- c("pseudotime", "celltype")
  pt_df <- pt_df[!is.na(pt_df$pseudotime), ]
  
  cat("=== ", pt_col, " by celltype ===\n", sep = "")
  stats <- pt_df %>% group_by(celltype) %>%
    summarise(n = n(),
              mean   = round(mean(pseudotime), 3),
              median = round(median(pseudotime), 3),
              sd     = round(sd(pseudotime), 3),
              min    = round(min(pseudotime), 3),
              max    = round(max(pseudotime), 3))
  print(stats)
  
  w <- wilcox.test(pseudotime ~ celltype, data = pt_df)
  cat("\nWilcoxon p-value:", signif(w$p.value, 4), "\n")
  
  p_box <- ggplot(pt_df, aes(x = celltype, y = pseudotime, fill = celltype)) +
    geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    stat_summary(fun = median, geom = "point", color = "red", size = 2) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = pt_col,
         title = paste0("CD14 vs CD16 (Wilcoxon p=", signif(w$p.value, 3), ")")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title  = element_text(hjust = 0.5, size = 10))
  
  p_dens <- ggplot(pt_df, aes(x = pseudotime, fill = celltype, color = celltype)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(x = pt_col, y = "Density", title = "Pseudotime distribution") +
    theme_classic() +
    theme(legend.position = "right",
          axis.text   = element_text(color = "black"),
          plot.title  = element_text(hjust = 0.5, size = 10))
  
  print(p_box)
  print(p_dens)
  invisible(list(stats = stats, wilcox = w, p_box = p_box, p_dens = p_dens))
}
# Slingshot
res_sling <- check_pseudotime(pbmc_mono, "sling_pseudotime")
# SCORPIUS
res_scorpius <- check_pseudotime(pbmc_mono, "scorpius_pseudotime")
# Monocle3
res_monocle3 <- check_pseudotime(pbmc_mono, "monocle3_pseudotime")
# DiffusionMap
res_diffusionmap<-check_pseudotime(pbmc_mono, "dpt_pseudotime")

#meta_sub<-pbmc_mono@meta.data
#saveRDS(meta_sub,file="E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\meta_sub.RDS")

######6 whole genome selection: chech the correlation between gene expression and pseudotime (at new server)
###### Load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)

###### Paths
work_dir <- "/mnt/data/work/masy/2026JIA"
input_rds <- file.path(work_dir, "pbmc_mono.RDS")
out_dir <- work_dir

###### Read data
pbmc_mono <- readRDS(input_rds)
cat("Cells:", ncol(pbmc_mono), "| Genes:", nrow(pbmc_mono), "\n")
cat("Pseudotime columns:",
    paste(grep("pseudotime", colnames(pbmc_mono@meta.data), value = TRUE), collapse = ", "), "\n")

###### Group annotation (3 groups)
pbmc_mono$group3 <- NA
pbmc_mono$group3[pbmc_mono$datasets %in% paste0("C", 1:6)] <- "HC"
pbmc_mono$group3[pbmc_mono$datasets %in% paste0("ploy", 1:6)] <- "Polygenic"
pbmc_mono$group3[pbmc_mono$datasets %in%
                   c("181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")] <- "Monogenic"
pbmc_mono$group3 <- factor(pbmc_mono$group3, levels = c("HC","Polygenic","Monogenic"))
cat("Groups:\n"); print(table(pbmc_mono$group3, useNA = "ifany"))
cat("Celltypes:\n"); print(table(pbmc_mono$celltype))

###### fit_gene_pt (memory-efficient)
fit_gene_pt <- function(seurat_obj, pt_col, gene_list, group_col = NULL,
                        cor_method = "spearman", min_cells = 30,
                        remove_zeros = FALSE, assay = "RNA", slot = "data") {
  if (!pt_col %in% colnames(seurat_obj@meta.data)) stop("pt_col not found")
  gene_list <- gene_list[gene_list %in% rownames(seurat_obj[[assay]])]
  if (length(gene_list) == 0) stop("No genes found")
  cat("Genes:", length(gene_list), "| Method:", cor_method,
      "| remove_zeros:", remove_zeros, "\n")
  pt <- seurat_obj@meta.data[[pt_col]]
  group <- if (!is.null(group_col)) seurat_obj@meta.data[[group_col]] else rep("all", ncol(seurat_obj))
  meta_df <- data.frame(pt = pt, group = group, row.names = colnames(seurat_obj))
  meta_df <- meta_df[!is.na(meta_df$pt), ]
  groups <- levels(droplevels(factor(meta_df$group)))
  cat("Groups:", paste(groups, collapse = ", "), "\n")
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  res_list <- list()
  for (g in groups) {
    cells_g <- rownames(meta_df)[meta_df$group == g]
    if (length(cells_g) < min_cells) { cat("Skip", g, ": n =", length(cells_g), "\n"); next }
    pt_g <- meta_df[cells_g, "pt"]
    for (gene in gene_list) {
      y <- as.numeric(expr_mat[gene, cells_g])
      x <- pt_g
      if (remove_zeros) { ok <- !is.na(x) & !is.na(y) & y > 0
      } else { ok <- !is.na(x) & !is.na(y) }
      if (sum(ok) < min_cells) next
      ct <- cor.test(x[ok], y[ok], method = cor_method)
      res_list[[length(res_list) + 1]] <- data.frame(
        gene = gene, group = g, n_cells = sum(ok),
        cor = as.numeric(ct$estimate), p_value = ct$p.value,
        stringsAsFactors = FALSE)
    }
  }
  res_df <- do.call(rbind, res_list) %>%
    group_by(group) %>% mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(group, p_adj) %>% as.data.frame()
  cat("\n=== Summary ===\n")
  print(res_df %>% group_by(group) %>%
          summarise(n_genes = n(), n_sig_p05 = sum(p_value < 0.05),
                    n_sig_adj05 = sum(p_adj < 0.05), .groups = "drop"))
  return(res_df)
}

###### run_one_pt_by_ct
run_one_pt_by_ct <- function(seurat_obj, pt_col, gene_list,
                             ct_col = "celltype",
                             ct_levels = c("CD14 Monocyte", "CD16 Monocyte"),
                             group_col = "group3",
                             cor_method = "spearman", min_cells = 30) {
  results <- list()
  for (ct in ct_levels) {
    ct_short <- gsub(" Monocyte", "", ct)
    cat("\n", strrep("-", 50), "\n")
    cat("Celltype:", ct, "\n")
    sub_obj <- seurat_obj[, seurat_obj@meta.data[[ct_col]] == ct]
    res_all <- fit_gene_pt(sub_obj, pt_col, gene_list, group_col = group_col,
                           cor_method = cor_method, min_cells = min_cells, remove_zeros = FALSE)
    res_all$version <- "all"; res_all$celltype <- ct_short
    res_no0 <- fit_gene_pt(sub_obj, pt_col, gene_list, group_col = group_col,
                           cor_method = cor_method, min_cells = min_cells, remove_zeros = TRUE)
    res_no0$version <- "no0"; res_no0$celltype <- ct_short
    results[[paste0(ct_short, "_all")]] <- res_all
    results[[paste0(ct_short, "_no0")]] <- res_no0
  }
  return(results)
}

###### filter_by_cor
filter_by_cor <- function(res_df, cor_threshold = 0.3) {
  genes_keep <- res_df %>% group_by(gene) %>%
    summarise(max_abs_cor = max(abs(cor), na.rm = TRUE), .groups = "drop") %>%
    filter(max_abs_cor > cor_threshold) %>% pull(gene)
  res_filtered <- res_df[res_df$gene %in% genes_keep, ]
  res_filtered <- res_filtered[order(res_filtered$gene, res_filtered$group), ]
  return(res_filtered)
}

###### output_excel
output_excel <- function(res_list, filename, cor_threshold = 0.3) {
  wb <- createWorkbook()
  for (nm in names(res_list)) {
    df <- filter_by_cor(res_list[[nm]], cor_threshold = cor_threshold)
    cat("Sheet", nm, ":", nrow(df), "rows (", length(unique(df$gene)), "genes )\n")
    addWorksheet(wb, nm)
    writeData(wb, nm, df)
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
  cat("Saved:", filename, "\n")
}

###### Genome-wide gene list
gene_list_all <- rownames(pbmc_mono)
cat("Total genes:", length(gene_list_all), "\n")

###### Run all 4 pseudotime methods
cat("\n", strrep("=", 60), "\n")
cat("Running Slingshot...\n")
res_sling <- run_one_pt_by_ct(pbmc_mono, "sling_pseudotime", gene_list_all)

cat("\n", strrep("=", 60), "\n")
cat("Running SCORPIUS...\n")
res_scorpius <- run_one_pt_by_ct(pbmc_mono, "scorpius_pseudotime", gene_list_all)

cat("\n", strrep("=", 60), "\n")
cat("Running Monocle3...\n")
res_monocle3 <- run_one_pt_by_ct(pbmc_mono, "monocle3_pseudotime", gene_list_all)

cat("\n", strrep("=", 60), "\n")
cat("Running DPT...\n")
res_dpt <- run_one_pt_by_ct(pbmc_mono, "dpt_pseudotime", gene_list_all)

###### Save full results as RDS
saveRDS(res_sling,    file.path(out_dir, "res_sling_genomewide.RDS"))
saveRDS(res_scorpius, file.path(out_dir, "res_scorpius_genomewide.RDS"))
saveRDS(res_monocle3, file.path(out_dir, "res_monocle3_genomewide.RDS"))
saveRDS(res_dpt,      file.path(out_dir, "res_dpt_genomewide.RDS"))
cat("\nFull RDS results saved.\n")

###### Output excel (|cor| > 0.3)
cat("\n", strrep("=", 60), "\n")
cat("Writing excel files...\n")
output_excel(res_sling,    file.path(out_dir, "gene_pt_sling_genomewide.xlsx"),    cor_threshold = 0.3)
output_excel(res_scorpius, file.path(out_dir, "gene_pt_scorpius_genomewide.xlsx"), cor_threshold = 0.3)
output_excel(res_monocle3, file.path(out_dir, "gene_pt_monocle3_genomewide.xlsx"), cor_threshold = 0.3)
output_excel(res_dpt,      file.path(out_dir, "gene_pt_dpt_genomewide.xlsx"),      cor_threshold = 0.3)

###### Quick summary
cat("\n", strrep("=", 60), "\n")
cat("=== Genes passing |cor| > 0.3 ===\n")
for (nm in c("CD14_all","CD14_no0","CD16_all","CD16_no0")) {
  cat(nm, ":\n")
  cat("  sling   :", length(unique(filter_by_cor(res_sling[[nm]])$gene)),    "genes\n")
  cat("  scorpius:", length(unique(filter_by_cor(res_scorpius[[nm]])$gene)), "genes\n")
  cat("  monocle3:", length(unique(filter_by_cor(res_monocle3[[nm]])$gene)), "genes\n")
  cat("  dpt     :", length(unique(filter_by_cor(res_dpt[[nm]])$gene)),      "genes\n")
}
cat("\nDone.\n")


#####7 loess-regression smooth
pbmc_mono<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\pbmc_mono.RDS")

###### Recreate group3 annotation
pbmc_mono$group3 <- NA
pbmc_mono$group3[pbmc_mono$datasets %in% paste0("C", 1:6)] <- "HC"
pbmc_mono$group3[pbmc_mono$datasets %in% paste0("ploy", 1:6)] <- "Polygenic"
pbmc_mono$group3[pbmc_mono$datasets %in%
                   c("181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")] <- "Monogenic"
pbmc_mono$group3 <- factor(pbmc_mono$group3, levels = c("HC","Polygenic","Monogenic"))
table(pbmc_mono$group3, useNA = "ifany")


###### filter_by_cor (extract genes from full RDS results)
filter_by_cor <- function(res_df, cor_threshold = 0.3) {
  library(dplyr)
  genes_keep <- res_df %>% group_by(gene) %>%
    summarise(max_abs_cor = max(abs(cor), na.rm = TRUE), .groups = "drop") %>%
    filter(max_abs_cor > cor_threshold) %>% pull(gene)
  res_filtered <- res_df[res_df$gene %in% genes_keep, ]
  res_filtered <- res_filtered[order(res_filtered$gene, res_filtered$group), ]
  return(res_filtered)
}

###### plot_gene_pt (loess curves for specified genes in one celltype)
plot_gene_pt <- function(seurat_obj, genes, pt_col, celltype,
                         version = "all",
                         colors = c("HC" = "#6699CC", "Polygenic" = "#FEC260",
                                    "Monogenic" = "#BB0021"),
                         group_col = "group3",
                         ct_col = "celltype",
                         assay = "RNA", slot = "data",
                         show_ci = TRUE, line_lwd = 0.8,
                         ncol = NULL) {
  library(ggplot2); library(Seurat); library(patchwork)
  if (!pt_col %in% colnames(seurat_obj@meta.data)) stop("pt_col '", pt_col, "' not found")
  if (!group_col %in% colnames(seurat_obj@meta.data)) stop("group_col '", group_col, "' not found")
  if (!ct_col %in% colnames(seurat_obj@meta.data)) stop("ct_col '", ct_col, "' not found")
  if (!celltype %in% unique(seurat_obj@meta.data[[ct_col]])) stop("celltype '", celltype, "' not found")
  sub_obj <- seurat_obj[, seurat_obj@meta.data[[ct_col]] == celltype]
  remove_zeros <- (version == "no0")
  cat("Celltype:", celltype, "| Version:", version,
      "| Genes:", length(genes), "| Cells:", ncol(sub_obj), "\n")
  p_list <- list()
  for (gene in genes) {
    if (!gene %in% rownames(sub_obj[[assay]])) { cat("Skip", gene, ": not found\n"); next }
    expr <- GetAssayData(sub_obj, assay = assay, slot = slot)[gene, ]
    df <- data.frame(
      pt = sub_obj@meta.data[[pt_col]],
      group = sub_obj@meta.data[[group_col]],
      expression = as.numeric(expr),
      row.names = colnames(sub_obj))
    df <- df[!is.na(df$pt) & !is.na(df$group) & !is.na(df$expression), ]
    df$group <- factor(df$group, levels = c("HC", "Polygenic", "Monogenic"))
    if (remove_zeros) df <- df[df$expression > 0, ]
    if (nrow(df) < 10) { cat("Skip", gene, ": too few cells\n"); next }
    title <- paste0(gene, " (", celltype, ", ", version, ")")
    p <- ggplot(df, aes(x = pt, y = expression, color = group, fill = group)) +
      geom_smooth(method = "loess", se = show_ci, linewidth = line_lwd,
                  alpha = ifelse(show_ci, 0.15, 0)) +
      coord_cartesian(ylim = c(0, NA)) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      labs(x = "Pseudotime", y = "Expression", title = title,
           color = NULL, fill = NULL) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, size = 9),
            legend.position = "right")
    p_list[[gene]] <- p
  }
  if (length(p_list) == 0) stop("No genes to plot")
  if (is.null(ncol)) ncol <- ceiling(sqrt(length(p_list)))
  combined <- wrap_plots(p_list, ncol = ncol)
  return(list(plots = p_list, combined = combined))
}


###### Read RDS, downstream
res_sling    <- readRDS("res_sling_genomewide.RDS")
res_scorpius <- readRDS("res_scorpius_genomewide.RDS")
res_monocle3 <- readRDS("res_monocle3_genomewide.RDS")
res_dpt      <- readRDS("res_dpt_genomewide.RDS")

###### Extract genes (|cor| > 0.5) from CD14_all of slingshot
genes_cd14_sling <- unique(filter_by_cor(res_sling$CD14_all, cor_threshold = 0.5)$gene)
cat("CD14_all slingshot genes:", length(genes_cd14_sling), "\n")

###### Slingshot (CD14_all,|cor| > 0.5)
p <- plot_gene_pt(pbmc_mono,
                  genes = genes_cd14_sling,
                  pt_col = "sling_pseudotime",
                  celltype = "CD14 Monocyte",
                  version = "all")
p$combined          # all genes combined
p$plots[[1]]        # first gene single plot


###### SCORPIUS (CD14_all, |cor| > 0.5)
genes_cd14_scorpius <- unique(filter_by_cor(res_scorpius$CD14_all, cor_threshold = 0.5)$gene)
cat("CD14_all scorpius genes:", length(genes_cd14_scorpius), "\n")
p_scorpius <- plot_gene_pt(pbmc_mono,
                           genes = genes_cd14_scorpius,
                           pt_col = "scorpius_pseudotime",
                           celltype = "CD14 Monocyte",
                           version = "all")
p_scorpius$combined

###### Monocle3 (CD14_all, |cor| > 0.5)
genes_cd14_monocle3 <- unique(filter_by_cor(res_monocle3$CD14_all, cor_threshold = 0.5)$gene)
cat("CD14_all monocle3 genes:", length(genes_cd14_monocle3), "\n")
p_monocle3 <- plot_gene_pt(pbmc_mono,
                           genes = genes_cd14_monocle3,
                           pt_col = "monocle3_pseudotime",
                           celltype = "CD14 Monocyte",
                           version = "all")
p_monocle3$combined

###### DPT (CD14_all, |cor| > 0.5)
genes_cd14_dpt <- unique(filter_by_cor(res_dpt$CD14_all, cor_threshold = 0.5)$gene)
cat("CD14_all dpt genes:", length(genes_cd14_dpt), "\n")
p_dpt <- plot_gene_pt(pbmc_mono,
                      genes = genes_cd14_dpt,
                      pt_col = "dpt_pseudotime",
                      celltype = "CD14 Monocyte",
                      version = "all")
p_dpt$combined


###### Slingshot (CD16_all, |cor| > 0.5)
genes_cd16_sling <- unique(filter_by_cor(res_sling$CD16_all, cor_threshold = 0.5)$gene)
cat("CD16_all slingshot genes:", length(genes_cd16_sling), "\n")
p_cd16_sling <- plot_gene_pt(pbmc_mono,
                             genes = genes_cd16_sling,
                             pt_col = "sling_pseudotime",
                             celltype = "CD16 Monocyte",
                             version = "all")
p_cd16_sling$combined

###### SCORPIUS (CD16_all, |cor| > 0.5)
genes_cd16_scorpius <- unique(filter_by_cor(res_scorpius$CD16_all, cor_threshold = 0.5)$gene)
cat("CD16_all scorpius genes:", length(genes_cd16_scorpius), "\n")
p_cd16_scorpius <- plot_gene_pt(pbmc_mono,
                                genes = genes_cd16_scorpius,
                                pt_col = "scorpius_pseudotime",
                                celltype = "CD16 Monocyte",
                                version = "all")
p_cd16_scorpius$combined

###### Monocle3 (CD16_all, |cor| > 0.5)
genes_cd16_monocle3 <- unique(filter_by_cor(res_monocle3$CD16_all, cor_threshold = 0.5)$gene)
cat("CD16_all monocle3 genes:", length(genes_cd16_monocle3), "\n")
p_cd16_monocle3 <- plot_gene_pt(pbmc_mono,
                                genes = genes_cd16_monocle3,
                                pt_col = "monocle3_pseudotime",
                                celltype = "CD16 Monocyte",
                                version = "all")
p_cd16_monocle3$combined

###### DPT (CD16_all, |cor| > 0.5)
genes_cd16_dpt <- unique(filter_by_cor(res_dpt$CD16_all, cor_threshold = 0.5)$gene)
cat("CD16_all dpt genes:", length(genes_cd16_dpt), "\n")
p_cd16_dpt <- plot_gene_pt(pbmc_mono,
                           genes = genes_cd16_dpt,
                           pt_col = "dpt_pseudotime",
                           celltype = "CD16 Monocyte",
                           version = "all")
p_cd16_dpt$combined


####print out pdf
###### output_plots_to_pdf: multi-page PDF (4x4 per page) from plot list
output_plots_to_pdf <- function(plot_list, output_file, ncol = 4, nrow = 4,
                                title = NULL, width = 16, height = 12) {
  library(patchwork)
  n_per_page <- ncol * nrow
  n_plots <- length(plot_list)
  n_pages <- ceiling(n_plots / n_per_page)
  cat("Total plots:", n_plots, "| Pages:", n_pages, "\n")
  pdf(output_file, width = width, height = height)
  for (i in seq_len(n_pages)) {
    start_idx <- (i - 1) * n_per_page + 1
    end_idx <- min(i * n_per_page, n_plots)
    page_plots <- plot_list[start_idx:end_idx]
    page_title <- if (!is.null(title)) paste0(title, " (page ", i, "/", n_pages, ")") else NULL
    p <- wrap_plots(page_plots, ncol = ncol) +
      plot_annotation(title = page_title,
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 12)))
    print(p)
  }
  dev.off()
  cat("Saved:", output_file, "\n")
}

###### Output all 8 results to PDF
output_plots_to_pdf(p$plots,              "CD14_sling_cor0.5.pdf",    title = "CD14 Monocyte - Slingshot (|cor| > 0.5)")
output_plots_to_pdf(p_scorpius$plots,     "CD14_scorpius_cor0.5.pdf", title = "CD14 Monocyte - SCORPIUS (|cor| > 0.5)")
output_plots_to_pdf(p_monocle3$plots,     "CD14_monocle3_cor0.5.pdf", title = "CD14 Monocyte - Monocle3 (|cor| > 0.5)")
output_plots_to_pdf(p_dpt$plots,           "CD14_dpt_cor0.5.pdf",      title = "CD14 Monocyte - DPT (|cor| > 0.5)")

output_plots_to_pdf(p_cd16_sling$plots,    "CD16_sling_cor0.5.pdf",    title = "CD16 Monocyte - Slingshot (|cor| > 0.5)")
output_plots_to_pdf(p_cd16_scorpius$plots, "CD16_scorpius_cor0.5.pdf", title = "CD16 Monocyte - SCORPIUS (|cor| > 0.5)")
output_plots_to_pdf(p_cd16_monocle3$plots, "CD16_monocle3_cor0.5.pdf", title = "CD16 Monocyte - Monocle3 (|cor| > 0.5)")
output_plots_to_pdf(p_cd16_dpt$plots,      "CD16_dpt_cor0.5.pdf",      title = "CD16 Monocyte - DPT (|cor| > 0.5)")

#######plot for scRNA figure
###### plot_gene_pt_multi (final: gene-only italic title, plain axis text)
plot_gene_pt_multi <- function(seurat_obj, gene_pt_map, celltype,
                               colors = c("HC" = "#6699CC", "Polygenic" = "#FEC260", "Monogenic" = "#AF478A"),
                               group_col = "group3",
                               ct_col = "celltype",
                               assay = "RNA", slot = "data",
                               show_ci = TRUE, line_lwd = 0.8, ci_alpha = 0.15) {
  library(ggplot2); library(Seurat)
  if (!group_col %in% colnames(seurat_obj@meta.data)) stop("group_col not found")
  if (!celltype %in% unique(seurat_obj@meta.data[[ct_col]])) stop("celltype not found")
  sub_obj <- seurat_obj[, seurat_obj@meta.data[[ct_col]] == celltype]
  pt_label_map <- c(sling_pseudotime = "Slingshot", scorpius_pseudotime = "SCORPIUS",
                    monocle3_pseudotime = "Monocle3", dpt_pseudotime = "DPT")
  cat("Celltype:", celltype, "| Genes:", length(gene_pt_map), "\n")
  result <- list()
  for (gene in names(gene_pt_map)) {
    if (!gene %in% rownames(sub_obj[[assay]])) { cat("Skip", gene, ": not found\n"); next }
    expr_full <- as.numeric(GetAssayData(sub_obj, assay = assay, slot = slot)[gene, ])
    pt_methods <- gene_pt_map[[gene]]
    gene_plots <- list()
    for (pt_col in pt_methods) {
      if (!pt_col %in% colnames(sub_obj@meta.data)) { cat("Skip", gene, pt_col, ": pt_col not found\n"); next }
      df <- data.frame(
        pt = sub_obj@meta.data[[pt_col]],
        group = as.character(sub_obj@meta.data[[group_col]]),
        expression = expr_full,
        stringsAsFactors = FALSE)
      df <- df[!is.na(df$pt) & !is.na(df$group) & !is.na(df$expression), ]
      df$group <- factor(df$group, levels = c("HC", "Polygenic", "Monogenic"))
      df <- df[!is.na(df$group), ]
      if (nrow(df) < 10) { cat("Skip", gene, pt_col, ": too few cells (", nrow(df), ")\n"); next }
      pt_label <- ifelse(pt_col %in% names(pt_label_map), pt_label_map[pt_col], pt_col)
      p <- ggplot(df, aes(x = pt, y = expression, color = group, fill = group)) +
        geom_smooth(method = "loess", se = show_ci, size = line_lwd,
                    alpha = ifelse(show_ci, ci_alpha, NA)) +
        coord_cartesian(ylim = c(0, NA)) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        labs(x = paste0(pt_label, " pseudotime"), y = "Expression",
             title = gene, color = NULL, fill = NULL) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              axis.line = element_blank(),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black"),
              plot.title = element_text(hjust = 0.5, size = 10, face = "italic"),
              legend.position = "none")
      gene_plots[[pt_col]] <- p
    }
    if (length(gene_plots) > 0) result[[gene]] <- gene_plots
  }
  cat("Generated plots for", length(result), "genes\n")
  return(result)
}

###### Main figure: multi-method consistent genes (3 methods each)
tier1_map <- list(
  "HIF1A" = c("dpt_pseudotime", "monocle3_pseudotime", "scorpius_pseudotime"),
  "TIMP1" = c("dpt_pseudotime", "monocle3_pseudotime", "scorpius_pseudotime")
)
plots_main <- plot_gene_pt_multi(pbmc_mono, tier1_map, celltype = "CD14 Monocyte")
HIF1A<-plots_main$HIF1A$dpt_pseudotime|plots_main$HIF1A$scorpius_pseudotime|plots_main$HIF1A$monocle3_pseudotime 
TIMP1<-plots_main$TIMP1$dpt_pseudotime|plots_main$TIMP1$scorpius_pseudotime|plots_main$TIMP1$monocle3_pseudotime
HIF1A/TIMP1
###### Supplementary figure: single-method specific genes
tier2_map <- list(
  "LGALS1"  = c("scorpius_pseudotime"),
  "SLC11A1" = c("scorpius_pseudotime"),
  "FTL"     = c("scorpius_pseudotime"),
  "TXNIP"   = c("dpt_pseudotime")
)
plots_supp <- plot_gene_pt_multi(pbmc_mono, tier2_map, celltype = "CD14 Monocyte")
sup<-plots_supp$LGALS1$scorpius_pseudotime | plots_supp$SLC11A1$scorpius_pseudotime | plots_supp$FTL$scorpius_pseudotime | plots_supp$TXNIP$dpt_pseudotime


################################################
################################################
########0916 revised  as single monogenic patient
#######1 loading data and group
###### Load packages
setwd("E:\\Cohort PPT\\JIA\\code\\CellTrajectory")
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
###### Read data
pbmc_mono <- readRDS("E:/Cohort PPT/JIA/code/CellTrajectory/pbmc_mono.RDS")
all_results <- readRDS("E:/Cohort PPT/JIA/code/CellTrajectory/genomewide_pt_mono_only.RDS")
###### Recreate group_mono annotation
pbmc_mono$group_mono <- NA
pbmc_mono$group_mono[pbmc_mono$datasets %in% paste0("C",1:6)] <- "HC"
mono_ids <- c("181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
pbmc_mono$group_mono[pbmc_mono$datasets %in% mono_ids] <- pbmc_mono$datasets[pbmc_mono$datasets %in% mono_ids]
pbmc_mono$group_mono <- factor(pbmc_mono$group_mono, levels = c("HC", mono_ids))
table(pbmc_mono$group_mono, useNA = "ifany")
###### filter_by_cor (adapted for new structure)
filter_by_cor <- function(res_df, cor_threshold = 0.4){
  genes_keep <- res_df %>% group_by(gene) %>%
    summarise(max_abs_cor = max(abs(cor), na.rm = TRUE), .groups = "drop") %>%
    filter(max_abs_cor > cor_threshold) %>% pull(gene)
  res_filtered <- res_df[res_df$gene %in% genes_keep, ]
  res_filtered <- res_filtered[order(res_filtered$gene, res_filtered$group), ]
  return(res_filtered)
}

#######2 ploting for 4 time and output
###### plot_gene_pt (batch: all filtered genes, one celltype, one pt method)
plot_gene_pt <- function(seurat_obj, genes, pt_col, celltype,
                         colors = c("HC"="#6699CC","181NOD2"="#CBDAA9","182LPIN2"="#E2A9C9",
                                    "183NOD2"="#B1DA99","190PSTPIP1"="#FFEFC1","200PSTPIP1"="#7B1FA2"),
                         group_col = "group_mono", ct_col = "celltype",
                         assay = "RNA", slot = "data",
                         show_ci = TRUE, line_lwd = 0.8, ncol = NULL){
  library(ggplot2); library(Seurat); library(patchwork)
  sub_obj <- seurat_obj[, seurat_obj@meta.data[[ct_col]] == celltype]
  cat("Celltype:", celltype, "| Genes:", length(genes), "| Cells:", ncol(sub_obj), "\n")
  p_list <- list()
  for(gene in genes){
    if(!gene %in% rownames(sub_obj[[assay]])){ cat("Skip", gene, "\n"); next }
    expr <- GetAssayData(sub_obj, assay = assay, slot = slot)[gene, ]
    df <- data.frame(pt = sub_obj@meta.data[[pt_col]],
                     group = as.character(sub_obj@meta.data[[group_col]]),
                     expression = as.numeric(expr), row.names = colnames(sub_obj))
    df <- df[!is.na(df$pt) & !is.na(df$group) & !is.na(df$expression), ]
    df$group <- factor(df$group, levels = names(colors))
    df <- df[!is.na(df$group), ]
    if(nrow(df) < 10) next
    p <- ggplot(df, aes(x = pt, y = expression, color = group, fill = group)) +
      geom_smooth(method = "loess", se = show_ci, linewidth = line_lwd,
                  alpha = ifelse(show_ci, 0.15, NA)) +
      coord_cartesian(ylim = c(0, NA)) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      labs(x = "Pseudotime", y = "Expression", title = gene, color = NULL, fill = NULL) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, size = 9, face = "italic"),
            legend.position = "right")
    p_list[[gene]] <- p
  }
  if(is.null(ncol)) ncol <- ceiling(sqrt(length(p_list)))
  combined <- wrap_plots(p_list, ncol = ncol)
  return(list(plots = p_list, combined = combined))
}

genes_cd14_sling <- unique(filter_by_cor(all_results$CD14$sling, 0.4)$gene)
cat("CD14 sling genes:", length(genes_cd14_sling), "\n")
p <- plot_gene_pt(pbmc_mono, genes = genes_cd14_sling,
                  pt_col = "sling_pseudotime", celltype = "CD14 Monocyte")
p$plots[[1]]

###### output_plots_to_pdf
output_plots_to_pdf <- function(plot_list, output_file, ncol = 4, nrow = 4,
                                title = NULL, width = 16, height = 12){
  n_per_page <- ncol * nrow
  n_plots <- length(plot_list)
  n_pages <- ceiling(n_plots / n_per_page)
  cat("Total plots:", n_plots, "| Pages:", n_pages, "\n")
  pdf(output_file, width = width, height = height)
  for(i in seq_len(n_pages)){
    start_idx <- (i-1) * n_per_page + 1
    end_idx <- min(i * n_per_page, n_plots)
    page_plots <- plot_list[start_idx:end_idx]
    p <- wrap_plots(page_plots, ncol = ncol) +
      plot_annotation(title = if(!is.null(title)) paste0(title, " (page ", i, "/", n_pages, ")") else NULL,
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 12)))
    print(p)
  }
  dev.off()
  cat("Saved:", output_file, "\n")
}

###### Batch: CD14
genes_cd14_sling <- unique(filter_by_cor(all_results$CD14$sling, 0.4)$gene)
p_cd14_sling <- plot_gene_pt(pbmc_mono, genes_cd14_sling, "sling_pseudotime", "CD14 Monocyte")
output_plots_to_pdf(p_cd14_sling$plots, "CD14_sling_cor04.pdf", title = "CD14 - Slingshot (|cor|>0.4)")

genes_cd14_scorpius <- unique(filter_by_cor(all_results$CD14$scorpius, 0.4)$gene)
p_cd14_scorpius <- plot_gene_pt(pbmc_mono, genes_cd14_scorpius, "scorpius_pseudotime", "CD14 Monocyte")
output_plots_to_pdf(p_cd14_scorpius$plots, "CD14_scorpius_cor04.pdf", title = "CD14 - SCORPIUS (|cor|>0.4)")

genes_cd14_monocle3 <- unique(filter_by_cor(all_results$CD14$monocle3, 0.4)$gene)
p_cd14_monocle3 <- plot_gene_pt(pbmc_mono, genes_cd14_monocle3, "monocle3_pseudotime", "CD14 Monocyte")
output_plots_to_pdf(p_cd14_monocle3$plots, "CD14_monocle3_cor04.pdf", title = "CD14 - Monocle3 (|cor|>0.4)")

genes_cd14_dpt <- unique(filter_by_cor(all_results$CD14$dpt, 0.4)$gene)
p_cd14_dpt <- plot_gene_pt(pbmc_mono, genes_cd14_dpt, "dpt_pseudotime", "CD14 Monocyte")
output_plots_to_pdf(p_cd14_dpt$plots, "CD14_dpt_cor04.pdf", title = "CD14 - DPT (|cor|>0.4)")

###### Batch: CD16
genes_cd16_sling <- unique(filter_by_cor(all_results$CD16$sling, 0.4)$gene)
p_cd16_sling <- plot_gene_pt(pbmc_mono, genes_cd16_sling, "sling_pseudotime", "CD16 Monocyte")
output_plots_to_pdf(p_cd16_sling$plots, "CD16_sling_cor04.pdf", title = "CD16 - Slingshot (|cor|>0.4)")

genes_cd16_scorpius <- unique(filter_by_cor(all_results$CD16$scorpius, 0.4)$gene)
p_cd16_scorpius <- plot_gene_pt(pbmc_mono, genes_cd16_scorpius, "scorpius_pseudotime", "CD16 Monocyte")
output_plots_to_pdf(p_cd16_scorpius$plots, "CD16_scorpius_cor04.pdf", title = "CD16 - SCORPIUS (|cor|>0.4)")

genes_cd16_monocle3 <- unique(filter_by_cor(all_results$CD16$monocle3, 0.4)$gene)
p_cd16_monocle3 <- plot_gene_pt(pbmc_mono, genes_cd16_monocle3, "monocle3_pseudotime", "CD16 Monocyte")
output_plots_to_pdf(p_cd16_monocle3$plots, "CD16_monocle3_cor04.pdf", title = "CD16 - Monocle3 (|cor|>0.4)")

genes_cd16_dpt <- unique(filter_by_cor(all_results$CD16$dpt, 0.4)$gene)
p_cd16_dpt <- plot_gene_pt(pbmc_mono, genes_cd16_dpt, "dpt_pseudotime", "CD16 Monocyte")
output_plots_to_pdf(p_cd16_dpt$plots, "CD16_dpt_cor04.pdf", title = "CD16 - DPT (|cor|>0.4)")


#####################################0922 adding plantir and for indepent monogenic patients
setwd("E:\\Cohort PPT\\JIA\\code\\CellTrajectory")
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
###### Read data
pbmc_mono <- readRDS("E:/Cohort PPT/JIA/code/CellTrajectory/pbmc_mono.RDS")
all_results <- readRDS("E:/Cohort PPT/JIA/code/CellTrajectory/genomewide_pt_mono_only.RDS")
###### Recreate group_mono annotation
pbmc_mono$group_mono <- NA
pbmc_mono$group_mono[pbmc_mono$datasets %in% paste0("C",1:6)] <- "HC"
mono_ids <- c("181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
pbmc_mono$group_mono[pbmc_mono$datasets %in% mono_ids] <- pbmc_mono$datasets[pbmc_mono$datasets %in% mono_ids]
pbmc_mono$group_mono <- factor(pbmc_mono$group_mono, levels = c("HC", mono_ids))
table(pbmc_mono$group_mono, useNA = "ifany")


###### filter_by_cor
filter_by_cor <- function(res_df, cor_threshold = 0.4){
  genes_keep <- res_df %>% group_by(gene) %>%
    summarise(max_abs_cor = max(abs(cor), na.rm = TRUE), .groups = "drop") %>%
    filter(max_abs_cor > cor_threshold) %>% pull(gene)
  res_df[res_df$gene %in% genes_keep, ]
}

###### plot_gene_pt
plot_gene_pt <- function(seurat_obj, genes, pt_col, celltype = NULL,
                         colors = c("HC"="#6699CC","181NOD2"="#CBDAA9","182LPIN2"="#E2A9C9",
                                    "183NOD2"="#B1DA99","190PSTPIP1"="#FFEFC1","200PSTPIP1"="#7B1FA2"),
                         group_col = "group_mono", ct_col = "celltype",
                         assay = "RNA", slot = "data",
                         show_ci = TRUE, line_lwd = 0.8, ncol = NULL){
  if(!is.null(celltype)){
    sub_obj <- seurat_obj[, seurat_obj@meta.data[[ct_col]] == celltype]
  } else {
    sub_obj <- seurat_obj
  }
  cat("Genes:", length(genes), "| Cells:", ncol(sub_obj), "\n")
  p_list <- list()
  for(gene in genes){
    if(!gene %in% rownames(sub_obj[[assay]])) next
    expr <- GetAssayData(sub_obj, assay = assay, slot = slot)[gene, ]
    df <- data.frame(pt = sub_obj@meta.data[[pt_col]],
                     group = as.character(sub_obj@meta.data[[group_col]]),
                     expression = as.numeric(expr), row.names = colnames(sub_obj))
    df <- df[!is.na(df$pt) & !is.na(df$group) & !is.na(df$expression), ]
    df$group <- factor(df$group, levels = names(colors))
    df <- df[!is.na(df$group), ]
    if(nrow(df) < 10) next
    p <- ggplot(df, aes(x = pt, y = expression, color = group, fill = group)) +
      geom_smooth(method = "loess", se = show_ci, linewidth = line_lwd,
                  alpha = ifelse(show_ci, 0.15, NA)) +
      coord_cartesian(ylim = c(0, NA)) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      labs(x = "Pseudotime", y = "Expression", title = gene, color = NULL, fill = NULL) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, size = 9, face = "italic"),
            legend.position = "right")
    p_list[[gene]] <- p
  }
  if(is.null(ncol)) ncol <- ceiling(sqrt(length(p_list)))
  combined <- wrap_plots(p_list, ncol = ncol)
  return(list(plots = p_list, combined = combined))
}

###### output_plots_to_pdf
output_plots_to_pdf <- function(plot_list, output_file, ncol = 4, nrow = 4,
                                title = NULL, width = 16, height = 12){
  n_per_page <- ncol * nrow
  n_plots <- length(plot_list)
  n_pages <- ceiling(n_plots / n_per_page)
  cat("Total plots:", n_plots, "| Pages:", n_pages, "\n")
  pdf(output_file, width = width, height = height)
  for(i in seq_len(n_pages)){
    start_idx <- (i-1) * n_per_page + 1
    end_idx <- min(i * n_per_page, n_plots)
    page_plots <- plot_list[start_idx:end_idx]
    p <- wrap_plots(page_plots, ncol = ncol) +
      plot_annotation(title = if(!is.null(title)) paste0(title, " (page ", i, "/", n_pages, ")") else NULL,
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 12)))
    print(p)
  }
  dev.off()
  cat("Saved:", output_file, "\n")
}

###### Check data structure
cat("CD14 methods:", paste(names(all_results$CD14), collapse = ", "), "\n")
cat("CD14_CD16 methods:", paste(names(all_results$CD14_CD16), collapse = ", "), "\n")

###### CD14
methods <- c("sling","scorpius","monocle3","dpt","palantir")
pt_cols <- c(sling="sling_pseudotime", scorpius="scorpius_pseudotime",
             monocle3="monocle3_pseudotime", dpt="dpt_pseudotime",
             palantir="palantir_pseudotime")
for(m in methods){
  if(!m %in% names(all_results$CD14)){ cat("Skip CD14", m, "\n"); next }
  genes <- unique(filter_by_cor(all_results$CD14[[m]], 0.4)$gene)
  cat("\nCD14", m, ":", length(genes), "genes\n")
  p <- plot_gene_pt(pbmc_mono, genes, pt_cols[m], celltype = "CD14 Monocyte")
  output_plots_to_pdf(p$plots, paste0("CD14_", m, "_cor04.pdf"),
                      title = paste0("CD14 - ", tools::toTitleCase(m), " (|cor|>0.4)"))
}

###### CD14_CD16 (no celltype filter)
for(m in methods){
  if(!m %in% names(all_results$CD14_CD16)){ cat("Skip CD14_CD16", m, "\n"); next }
  genes <- unique(filter_by_cor(all_results$CD14_CD16[[m]], 0.4)$gene)
  cat("\nCD14_CD16", m, ":", length(genes), "genes\n")
  p <- plot_gene_pt(pbmc_mono, genes, pt_cols[m], celltype = NULL)
  output_plots_to_pdf(p$plots, paste0("CD14CD16_", m, "_cor04.pdf"),
                      title = paste0("CD14+CD16 - ", tools::toTitleCase(m), " (|cor|>0.4)"))
}

######check data
meta <- pbmc_mono@meta.data
meta <- meta[meta$celltype %in% c("CD14 Monocyte","CD16 Monocyte"), ]
pt_methods <- c("sling_pseudotime","scorpius_pseudotime","monocle3_pseudotime","dpt_pseudotime","palantir_pseudotime")
patients <- c("181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
res_list <- list()
for(pt in pt_methods){
  for(ptnt in patients){
    sub <- meta[meta$datasets == ptnt, ]
    if(nrow(sub) < 10) next
    med14 <- median(sub[[pt]][sub$celltype == "CD14 Monocyte"], na.rm = TRUE)
    med16 <- median(sub[[pt]][sub$celltype == "CD16 Monocyte"], na.rm = TRUE)
    res_list[[length(res_list)+1]] <- data.frame(
      method = gsub("_pseudotime","",pt),
      patient = ptnt,
      median_CD14 = round(med14, 4),
      median_CD16 = round(med16, 4),
      CD16_higher = med16 > med14
    )
  }
}
dir_df <- do.call(rbind, res_list)
dir_df<-dir_df[dir_df$CD16_higher=="TRUE",]


###### Plot gene across valid patient×method combinations
###### Batch: boxplot + loess for all valid combos
plot_gene_batch <- function(seurat_obj, gene, valid_df,
                            ct_colors = c("CD14 Monocyte"="#6699CC","CD16 Monocyte"="#F49D5C"),
                            line_color = "#6699CC",
                            assay = "RNA", slot = "data"){
  library(ggplot2); library(patchwork); library(dplyr)
  pt_col_map <- c(sling="sling_pseudotime", scorpius="scorpius_pseudotime",
                  monocle3="monocle3_pseudotime", dpt="dpt_pseudotime",
                  palantir="palantir_pseudotime")
  valid_df <- valid_df[valid_df$CD16_higher == TRUE, ]
  box_list <- list()
  loess_list <- list()
  for(i in seq_len(nrow(valid_df))){
    ptm <- valid_df$method[i]
    ptnt <- valid_df$patient[i]
    pt_col <- pt_col_map[ptm]
    sub <- seurat_obj[, seurat_obj$datasets == ptnt]
    meta <- sub@meta.data
    meta <- meta[!is.na(meta[[pt_col]]), ]
    meta$celltype <- factor(meta$celltype, levels = names(ct_colors))
    meta <- meta[!is.na(meta$celltype), ]
    p_box <- ggplot(meta, aes(x = celltype, y = .data[[pt_col]], fill = celltype)) +
      geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +
      scale_fill_manual(values = ct_colors) +
      labs(title = paste0(ptnt,"|",ptm), x = "", y = "") +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 7),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
            axis.text.y = element_text(color = "black", size = 6))
    box_list[[paste0(ptnt,"_",ptm)]] <- p_box
    expr <- GetAssayData(sub, assay = assay, slot = slot)[gene, ]
    df2 <- data.frame(pt = sub@meta.data[[pt_col]],
                      expression = as.numeric(expr))
    df2 <- df2[!is.na(df2$pt) & !is.na(df2$expression), ]
    p_loess <- ggplot(df2, aes(x = pt, y = expression)) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, alpha = 0.15, color = line_color, fill = line_color) +
      labs(title = paste0(ptnt,"|",ptm), x = "", y = "") +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 7),
            axis.text = element_text(color = "black", size = 6))
    loess_list[[paste0(ptnt,"_",ptm)]] <- p_loess
  }
  ncol <- ceiling(sqrt(length(box_list)))
  out <- list(
    boxplot = wrap_plots(box_list, ncol = ncol),
    loess = wrap_plots(loess_list, ncol = ncol)
  )
  return(out)
}

res <- plot_gene_batch(pbmc_mono, IFN1$gene[33], dir_df[dir_df$patient=="181NOD2",])
res$loess 
#res$boxplot
c("TRIM16","DDX60","IFIT5")#182LPIN2




list(
  c("182LPIN2","IFI30","SCROPIUS") ,
  c("200PSTPIP1","IFI30","SCROPIUS"),
  c("182LPIN2","LGALS3","SCROPIUS") ,
  c("200PSTPIP1","LGALS3","SCROPIUS"),
  c("182LPIN2","HLA-DRA","SCROPIUS") ,
  c("200PSTPIP1","HLA-DRA","SCROPIUS"),
)


###### Batch: loess curves per patient per gene (20 per page)
plot_gene_pdf <- function(seurat_obj, gene_vec, valid_df, out_dir = ".",
                          line_color = "#6699CC",
                          assay = "RNA", slot = "data"){
  library(ggplot2); library(dplyr); library(patchwork)
  pt_col_map <- c(sling="sling_pseudotime", scorpius="scorpius_pseudotime",
                  monocle3="monocle3_pseudotime", dpt="dpt_pseudotime",
                  palantir="palantir_pseudotime")
  valid_df <- valid_df[valid_df$CD16_higher == TRUE, ]
  patients <- unique(valid_df$patient)
  for(ptnt in patients){
    sub_df <- valid_df[valid_df$patient == ptnt, ]
    methods <- sub_df$method
    methods <- methods[methods %in% names(pt_col_map)]
    if(length(methods) == 0) next
    sub <- seurat_obj[, seurat_obj$datasets == ptnt]
    p_list <- list()
    idx <- 0
    for(gene in gene_vec){
      if(!gene %in% rownames(sub[[assay]])) next
      expr <- GetAssayData(sub, assay = assay, slot = slot)[gene, ]
      for(ptm in methods){
        pt_col <- pt_col_map[ptm]
        df2 <- data.frame(pt = sub@meta.data[[pt_col]],
                          expression = as.numeric(expr))
        df2 <- df2[!is.na(df2$pt) & !is.na(df2$expression), ]
        p <- ggplot(df2, aes(x = pt, y = expression)) +
          geom_smooth(method = "loess", se = TRUE, linewidth = 0.7, alpha = 0.15,
                      color = line_color, fill = line_color) +
          labs(title = paste0(gene, " | ", ptm), x = "", y = "") +
          theme_classic() +
          theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
                axis.line = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 8, face = "italic"),
                axis.text = element_text(color = "black", size = 6))
        idx <- idx + 1
        p_list[[idx]] <- p
      }
    }
    if(length(p_list) == 0) next
    n_per_page <- 20
    n_pages <- ceiling(length(p_list) / n_per_page)
    pdf(file.path(out_dir, paste0(ptnt, "_genes_loess.pdf")), width = 16, height = 18)
    for(pg in seq_len(n_pages)){
      start_idx <- (pg-1)*n_per_page + 1
      end_idx <- min(pg*n_per_page, length(p_list))
      page_plots <- p_list[start_idx:end_idx]
      print(wrap_plots(page_plots, ncol = 4))
    }
    dev.off()
    cat("Saved:", ptnt, "_genes_loess.pdf", "|", length(p_list), "plots |", n_pages, "pages\n")
  }
}

plot_gene_pdf(pbmc_mono, gene_vec = c("LYZ","S100A8","S100A9","FCGR3B","CD14","VCAN",IFN1$gene),
              valid_df = dir_df, out_dir = ".")

#####plot for figures
###### Final: boxplot + gene loess per timepoint
###### Final: boxplot + gene loess per timepoint
plot_pt_gene_final <- function(seurat_obj, patient, gene_vec, pt_cols,
                               curve_color = "#E2A9C9",
                               ct_colors = c("CD14 Monocyte"="#AF478A","CD16 Monocyte"="#EC706E"),
                               assay = "RNA", slot = "data"){
  library(ggplot2); library(patchwork)
  sub <- seurat_obj[, seurat_obj$datasets == patient]
  cat("Patient:", patient, "| Cells:", ncol(sub), "\n")
  row_list <- list()
  n_rows <- length(pt_cols)
  for(r in seq_along(pt_cols)){
    pt_col <- pt_cols[r]
    meta <- sub@meta.data
    meta <- meta[!is.na(meta[[pt_col]]), ]
    meta$celltype <- factor(meta$celltype, levels = names(ct_colors))
    meta <- meta[!is.na(meta$celltype), ]
    wt <- wilcox.test(meta[[pt_col]][meta$celltype == "CD14 Monocyte"],
                      meta[[pt_col]][meta$celltype == "CD16 Monocyte"])
    p_str <- ifelse(wt$p.value < 2.2e-16, "<2.2e-16", formatC(wt$p.value, format = "e", digits = 2))
    vals <- meta[[pt_col]]
    ymin <- min(vals, na.rm = TRUE); ymax <- max(vals, na.rm = TRUE)
    y_range <- ymax - ymin
    y_bracket <- ymax + 0.08*y_range
    is_last <- (r == n_rows)
    p_box <- ggplot(meta, aes(x = celltype, y = .data[[pt_col]], fill = celltype)) +
      geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +
      scale_fill_manual(values = ct_colors) +
      scale_y_continuous(limits = c(ymin - 0.02*y_range, y_bracket + 0.12*y_range), expand = c(0,0)) +
      annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, color = "black", linewidth = 0.5) +
      annotate("segment", x = 1, xend = 1, y = y_bracket - 0.02*y_range, yend = y_bracket, color = "black", linewidth = 0.5) +
      annotate("segment", x = 2, xend = 2, y = y_bracket - 0.02*y_range, yend = y_bracket, color = "black", linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_bracket + 0.02*y_range, label = paste0("P=",p_str), size = 3, fontface = "italic") +
      labs(title = pt_col, x = "", y = "") +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 8),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
            axis.ticks.x = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size = 7),
            plot.margin = margin(2, 2, 2, 2))
    if(!is_last){
      p_box <- p_box + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    gene_plots <- list()
    for(gene in gene_vec){
      if(!gene %in% rownames(sub[[assay]])) next
      expr <- GetAssayData(sub, assay = assay, slot = slot)[gene, ]
      df2 <- data.frame(pt = sub@meta.data[[pt_col]],
                        expression = as.numeric(expr))
      df2 <- df2[!is.na(df2$pt) & !is.na(df2$expression), ]
      p_g <- ggplot(df2, aes(x = pt, y = expression)) +
        geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, alpha = 0.15,
                    color = curve_color, fill = curve_color) +
        labs(title = gene, x = "", y = "") +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
              axis.line = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 8, face = "italic"),
              axis.text = element_text(color = "black", size = 7),
              plot.margin = margin(2, 2, 2, 2))
      if(!is_last){
        p_g <- p_g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      gene_plots[[gene]] <- p_g
    }
    row_combined <- p_box | wrap_plots(gene_plots, ncol = length(gene_plots))
    row_combined <- row_combined + plot_layout(widths = c(0.6, length(gene_plots)))
    row_list[[pt_col]] <- row_combined
  }
  final <- wrap_plots(row_list, ncol = 1) +
    plot_layout(heights = rep(1, n_rows))
  return(final)
}

p <- plot_pt_gene_final(pbmc_mono, patient = "182LPIN2",
                        gene_vec = c("RSAD2","IFIT5","DDX58"),
                        pt_cols = c("sling_pseudotime","scorpius_pseudotime"),
                        curve_color = "#E2A9C9")
p
