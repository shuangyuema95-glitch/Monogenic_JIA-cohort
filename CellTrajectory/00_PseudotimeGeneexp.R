library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(SCORPIUS)
library(ggplot2)
library(dplyr)


setwd("E:\\Cohort PPT\\JIA\\code\\CellTrajectory")
# pbmc1<-readRDS("E:\\Cohort PPT\\JIA\\code\\pbmc1.RDS")
# pbmc_mono<-subset(pbmc1,celltype%in%c("CD14 Monocyte","CD16 Monocyte"))
pbmc_mono<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\pbmc_mono.RDS")

##########Step1 pseudotime calculation and fit regression
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
length(pbmc_mono$sling_pseudotime)#25857
ncol(pbmc_mono)#25857
#rm(sce, sce_sling, sling_pt); gc()

######2 SCORPIUS
library(SCORPIUS)
expr_mat <- as.matrix(t(GetAssayData(pbmc_mono, assay = "RNA", slot = "data")))
space <- reduce_dimensionality(expr_mat, "spearman")
traj  <- infer_trajectory(space)
pbmc_mono$scorpius_pseudotime <- traj$time[colnames(pbmc_mono)]
sum(is.na(pbmc_mono$scorpius_pseudotime))#0
length(pbmc_mono$scorpius_pseudotime)#25857
ncol(pbmc_mono)#25857

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

# ######4 DiffusionMap
# dpt<-read.csv("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\dpt_pseudotime.csv")
# ncol(pbmc_mono)
# dim(dpt)
# pbmc_mono@meta.data$dpt_pseudotime<-dpt[match(rownames(pbmc_mono@meta.data),dpt$X),"dpt_pseudotime"]
# 
# 
# ######5 plantir
# pal<-read.csv("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\palantir_pseudotime.csv")
# ncol(pbmc_mono)
# dim(pal)
# pbmc_mono@meta.data$palantir_pseudotime<-pal[match(rownames(pbmc_mono@meta.data),pal$X),"X0"]
# sum(is.na(pbmc_mono$palantir_pseudotime))


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

res_sling <- check_pseudotime(pbmc_mono, "sling_pseudotime")# Slingshot
res_scorpius <- check_pseudotime(pbmc_mono, "scorpius_pseudotime")# SCORPIUS
res_monocle3 <- check_pseudotime(pbmc_mono, "monocle3_pseudotime")# Monocle3
#meta_sub<-pbmc_mono@meta.data
#saveRDS(meta_sub,file="E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\meta_sub.RDS")

######6 whole genome selection: chech the correlation between gene expression and pseudotime (at new server)
###### Load packages
library(Seurat); library(dplyr); library(tidyr); library(openxlsx)
###### Paths
work_dir <- "/mnt/data/work/masy/2026JIA"
input_rds <- file.path(work_dir, "pbmc_mono.RDS")
out_dir <- work_dir
###### Read data
pbmc_mono <- readRDS(input_rds)
cat("Cells:", ncol(pbmc_mono), "| Genes:", nrow(pbmc_mono), "\n")
cat("Pseudotime columns:",
    paste(grep("pseudotime", colnames(pbmc_mono@meta.data), value = TRUE), collapse = ", "), "\n")
###### Group annotation (HC + 6 monogenic separately)
pbmc_mono$group_mono <- NA
pbmc_mono$group_mono[pbmc_mono$datasets %in% paste0("C",1:6)] <- "HC"
mono_ids <- c("157LACC1","181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
pbmc_mono$group_mono[pbmc_mono$datasets %in% mono_ids] <- pbmc_mono$datasets[pbmc_mono$datasets %in% mono_ids]
pbmc_mono$group_mono <- factor(pbmc_mono$group_mono, levels = c("HC", mono_ids))
cat("Groups:\n"); print(table(pbmc_mono$group_mono, useNA = "ifany"))
cat("Celltypes:\n"); print(table(pbmc_mono$celltype))
###### fit_gene_pt (no remove_zeros)
fit_gene_pt <- function(seurat_obj, pt_col, gene_list, group_col = "group_mono",
                        cor_method = "spearman", min_cells = 30,
                        assay = "RNA", slot = "data"){
  if(!pt_col %in% colnames(seurat_obj@meta.data)) stop("pt_col not found")
  gene_list <- gene_list[gene_list %in% rownames(seurat_obj[[assay]])]
  if(length(gene_list) == 0) stop("No genes found")
  cat("Genes:", length(gene_list), "| Method:", cor_method, "\n")
  pt <- seurat_obj@meta.data[[pt_col]]
  group <- seurat_obj@meta.data[[group_col]]
  meta_df <- data.frame(pt = pt, group = group, row.names = colnames(seurat_obj))
  meta_df <- meta_df[!is.na(meta_df$pt) & !is.na(meta_df$group), ]
  groups <- levels(droplevels(factor(meta_df$group)))
  cat("Groups:", paste(groups, collapse = ", "), "\n")
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  res_list <- list()
  for(g in groups){
    cells_g <- rownames(meta_df)[meta_df$group == g]
    if(length(cells_g) < min_cells){ cat("Skip", g, ": n =", length(cells_g), "\n"); next }
    pt_g <- meta_df[cells_g, "pt"]
    for(gene in gene_list){
      y <- as.numeric(expr_mat[gene, cells_g])
      x <- pt_g
      ok <- !is.na(x) & !is.na(y)
      if(sum(ok) < min_cells) next
      ct <- cor.test(x[ok], y[ok], method = cor_method)
      res_list[[length(res_list)+1]] <- data.frame(
        gene = gene, group = g, n_cells = sum(ok),
        cor = as.numeric(ct$estimate), p_value = ct$p.value, stringsAsFactors = FALSE)
    }
  }
  res_df <- do.call(rbind, res_list) %>%
    group_by(group) %>% mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(group, p_adj) %>% as.data.frame()
  cat("Summary:\n")
  print(res_df %>% group_by(group) %>%
          summarise(n_genes = n(), n_sig = sum(p_adj < 0.05), .groups = "drop"))
  return(res_df)
}
###### run_subset (one subset, filter low expr, pt methods)
run_subset <- function(seurat_obj, sub_name, sub_cells, pt_methods, min_expr_cells = 10,
                       min_cells = 30){
  cat("\n", strrep("=", 60), "\n")
  cat("Subset:", sub_name, "\n")
  sub_obj <- seurat_obj[, sub_cells]
  expr_mat <- GetAssayData(sub_obj, assay = "RNA", slot = "data")
  n_expr <- rowSums(expr_mat > 0)
  genes_keep <- names(n_expr)[n_expr >= min_expr_cells]
  cat("Genes total:", length(n_expr), "| expressed in >=", min_expr_cells, "cells:", length(genes_keep), "\n")
  results <- list()
  for(ptm in pt_methods){
    pt_short <- gsub("_pseudotime", "", ptm)
    cat("\n--- ", ptm, " ---\n")
    res <- fit_gene_pt(sub_obj, ptm, genes_keep, min_cells = min_cells)
    res$subset <- sub_name
    results[[pt_short]] <- res
  }
  return(results)
}
###### Run
pt_methods <- c("sling_pseudotime","scorpius_pseudotime","monocle3_pseudotime","dpt_pseudotime","palantir_pseudotime")
pt_methods <- intersect(pt_methods, colnames(pbmc_mono@meta.data))
cat("Methods to run:", paste(pt_methods, collapse = ", "), "\n")
all_results <- list()
all_results[["CD14"]] <- run_subset(pbmc_mono, "CD14",
                                    pbmc_mono$celltype == "CD14 Monocyte",
                                    pt_methods, min_expr_cells = 10, min_cells = 30)
all_results[["CD14_CD16"]] <- run_subset(pbmc_mono, "CD14_CD16",
                                         rep(TRUE, ncol(pbmc_mono)),
                                         pt_methods, min_expr_cells = 10, min_cells = 30)
###### Save RDS
saveRDS(all_results, file.path(out_dir, "genomewide_pt_mono_only.RDS"))
cat("\nSaved: genomewide_pt_mono_only.RDS\n")
###### Output excel (|cor| > 0.4)
wb <- createWorkbook()
for(sub in names(all_results)){
  for(pt in names(all_results[[sub]])){
    df <- all_results[[sub]][[pt]]
    df <- df[abs(df$cor) > 0.4, ]
    sheet_name <- paste0(sub, "_", pt)
    cat("Sheet", sheet_name, ":", nrow(df), "rows (", length(unique(df$gene)), "genes)\n")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
  }
}
saveWorkbook(wb, file.path(out_dir, "gene_pt_mono_only_cor04.xlsx"), overwrite = TRUE)
cat("Saved: gene_pt_mono_only_cor04.xlsx\n")
###### Summary
cat("\n=== Genes passing |cor| > 0.4 ===\n")
for(sub in names(all_results)){
  cat(sub, ":\n")
  for(pt in names(all_results[[sub]])){
    df <- all_results[[sub]][[pt]]
    cat("  ", pt, ":", length(unique(df$gene[abs(df$cor) > 0.4])), "genes\n")
  }
}
cat("\nDone.\n")

##########Step2 fit plot for donwstream
#####################################0922 adding plantir and for indepent monogenic patients
setwd("E:\\Cohort PPT\\JIA\\code\\CellTrajectory")
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
######1 Read data
pbmc_mono <- readRDS("E:/Cohort PPT/JIA/code/CellTrajectory/pbmc_mono.RDS")
all_results <- readRDS("E:/Cohort PPT/JIA/code/CellTrajectory/genomewide_pt_mono_only.RDS")
table(pbmc_mono$datasets) # 923 polygenic
######2 Recreate group_mono annotation
pbmc_mono$group_mono <- NA
pbmc_mono$group_mono[pbmc_mono$datasets %in% paste0("C",1:6)] <- "HC"
mono_ids <- c("157LACC1","181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
pbmc_mono$group_mono[pbmc_mono$datasets %in% mono_ids] <- pbmc_mono$datasets[pbmc_mono$datasets %in% mono_ids]
pbmc_mono$group_mono <- factor(pbmc_mono$group_mono, levels = c("HC", mono_ids))
table(pbmc_mono$group_mono, useNA = "ifany") #93 NA

######3 check data
meta <- pbmc_mono@meta.data
meta <- meta[meta$celltype %in% c("CD14 Monocyte","CD16 Monocyte"), ]
pt_methods <- c("sling_pseudotime","scorpius_pseudotime","monocle3_pseudotime")
patients <- c("157LACC1","181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
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


######4 Plot gene across valid patient×method combinations
####(1) Batch: boxplot + loess for all valid combos
IFN1 <- c("DDX60","EPSTI1","HERC5","HERC6","IFI27","IFI44","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFIT5","ISG15","LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SIGLEC1","SPATS2L","USP18","CXCL10","GBP1","SOCS1")
IFN2<- c("GPR18","TRIM16","CXCL9","SLCSA6","GPR146","CIITA","PRF1","JAK1","IFNG","ATP23","B4GALT3","TRIB3","CCR5","DNA2","AHR","IDO1","SELP","JAK2","RAPGEF6","BMF","STAT1","FBXL5","IRF1","TAP1","IFITM1","GBP2","IFIH1")


plot_gene_batch <- function(seurat_obj, gene, valid_df,
                            ct_colors = c("CD14 Monocyte"="#6699CC","CD16 Monocyte"="#F49D5C"),
                            line_color = "#6699CC",
                            assay = "RNA", slot = "data"){
  library(ggplot2); library(patchwork); library(dplyr)
  pt_col_map <- c(sling="sling_pseudotime", scorpius="scorpius_pseudotime",
                  monocle3="monocle3_pseudotime")
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

res <- plot_gene_batch(pbmc_mono, IFN1[1], dir_df[dir_df$patient=="181NOD2",])
res$loess 
res$boxplot


####(2) Batch: loess curves per patient per gene (20 per page)
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

plot_gene_pdf(pbmc_mono, gene_vec = c("LYZ","S100A8","S100A9","FCGR3B","CD14","VCAN",IFN1),
              valid_df = dir_df, out_dir = ".")

######5 Plot gene across valid patient×method combinations
#####(1) Final: boxplot + gene loess per timepoint
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
      gene_plots[[gene]] <- p_g
    }
    row_combined <- p_box | wrap_plots(gene_plots, ncol = length(gene_plots))
    row_combined <- row_combined + plot_layout(widths = c(0.5, length(gene_plots)))
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

#####(2) increase tread for each sample
PT_COL_MAP <- c(sling = "sling_pseudotime", scorpius = "scorpius_pseudotime", monocle3 = "monocle3_pseudotime")
PATIENT_COLORS <- c("HC" = "#6699CC", "157LACC1" = "#7E75AB", "181NOD2" = "#CBDAA9", "182LPIN2" = "#E2A9C9",
                    "183NOD2" = "#B1DA99", "190PSTPIP1" = "#FFEFC1", "200PSTPIP1" = "#7B1FA2")
GENE_PAL <- c("#0B71AB","#EC706E","#7E75AB","#86C7B4","#FEC260","#AF478A","#6699CC","#F49D5C",
              "#847AB3","#7B1FA2","#4DAF4A","#9CD2ED","#D9C2D9","#F1BAAB","#6BA292","#CBD281")

pt_gene_hits <- list(
  "182LPIN2"   = list(scorpius = c("EPSTI1","IFI27","IFI44","IFI44L","IFIT1","IFIT2","LY6E")),
  "181NOD2"    = list(sling = "DDX60"),
  "157LACC1"   = list(scorpius = c("IFI27","IFI44"), sling = "OAS1"),
  "183NOD2"    = list(monocle3 = c("IFI44","IFIT3","GBP1"), sling = "IFIT5"),
  "200PSTPIP1" = list(sling = c("DDX60","IFIT6","IFI44","IFI44L")),
  "190PSTPIP1" = list(sling = c("CXCL10","LY6E","OAS1","IFI6","IFIT1","EPSTI1"))
)

pt_gene_df <- function(seurat_obj, hits, patient, assay = "RNA", slot = "data",
                       pt_map = PT_COL_MAP, normalize = FALSE, min_cells = 10){
  sub <- seurat_obj[, seurat_obj$datasets %in% patient]
  if(ncol(sub) == 0){ warning("no cells: ", patient); return(NULL) }
  md <- sub@meta.data
  out <- vector("list", 0)
  k <- 0
  for(m in names(hits[[patient]])){
    pt_col <- unname(pt_map[m])
    if(is.na(pt_col) || !pt_col %in% colnames(md)){ warning("bad method: ", patient, "/", m); next }
    genes <- hits[[patient]][[m]]
    miss <- setdiff(genes, rownames(sub[[assay]]))
    if(length(miss)) warning(patient, "/", m, " absent in assay: ", paste(miss, collapse = ","))
    genes <- intersect(genes, rownames(sub[[assay]]))
    if(!length(genes)) next
    mat <- GetAssayData(sub, assay = assay, slot = slot)[genes, , drop = FALSE]
    pt_all <- md[[pt_col]]
    for(g in genes){
      v <- as.numeric(mat[g, ])
      ok <- !is.na(pt_all) & !is.na(v)
      if(sum(ok) < min_cells){ warning("too few cells: ", patient, "/", g); next }
      x <- pt_all[ok]
      if(normalize){
        r <- range(x)
        if(diff(r) == 0) next
        x <- (x - r[1]) / diff(r)
      }
      k <- k + 1
      out[[k]] <- data.frame(patient = patient, method = m, gene = g, pt = x, expression = v[ok])
    }
  }
  if(!k) return(NULL)
  do.call(rbind, out)
}

plot_pt_gene_by_patient <- function(seurat_obj, hits, patients = names(hits), ncol = 3, nrow = NULL,
                                    ci = TRUE, pt_map = PT_COL_MAP, pt_cols = PATIENT_COLORS,
                                    assay = "RNA", slot = "data",
                                    xlab = "Pseudotime", ylab = "Expression"){
  library(ggplot2)
  res <- list()
  for(p in patients){
    df <- pt_gene_df(seurat_obj, hits, p, assay, slot, pt_map, normalize = FALSE)
    if(is.null(df)) next
    multi <- length(unique(df$method)) > 1
    lab <- if(multi) paste0(df$gene, " (", df$method, ")") else df$gene
    df$panel <- factor(lab, levels = unique(lab))
    col <- if(p %in% names(pt_cols)) unname(pt_cols[p]) else "#6699CC"
    res[[p]] <- ggplot(df, aes(x = pt, y = expression)) +
      geom_smooth(method = "loess", se = ci, linewidth = 0.8, alpha = 0.15, color = col, fill = col) +
      facet_wrap(~panel, ncol = ncol, nrow = nrow, scales = "free") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, NA)) +
      labs(title = p, x = xlab, y = ylab) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9),
            strip.background = element_blank(),
            strip.text = element_text(size = 7, face = "italic", color = "black"),
            axis.text = element_text(color = "black", size = 7))
  }
  res
}

plot_pt_gene_overlay <- function(seurat_obj, hits, patients = names(hits), gene_colors = NULL,
                                 normalize = TRUE, span = 0.75, label_size = 2.2,
                                 pt_map = PT_COL_MAP, assay = "RNA", slot = "data",
                                 xlab = "Pseudotime (scaled 0-1)", ylab = "Expression"){
  library(ggplot2)
  if(is.null(gene_colors)){
    gs <- unique(unlist(lapply(hits[patients], function(x) unlist(x, use.names = FALSE))))
    gene_colors <- setNames(rep(GENE_PAL, length.out = length(gs)), gs)
  }
  res <- list()
  for(p in patients){
    df <- pt_gene_df(seurat_obj, hits, p, assay, slot, pt_map, normalize = normalize)
    if(is.null(df)) next
    multi <- length(unique(df$method)) > 1
    df$key <- if(multi) paste0(df$gene, " (", df$method, ")") else df$gene
    fit <- do.call(rbind, lapply(split(df, df$key), function(d){
      if(nrow(d) < 10 || length(unique(d$pt)) < 5) return(NULL)
      f <- loess(expression ~ pt, data = d, span = span)
      gx <- seq(min(d$pt), max(d$pt), length.out = 200)
      data.frame(key = d$key[1], gene = d$gene[1], pt = gx, fit = predict(f, gx))
    }))
    if(is.null(fit)){ warning("loess failed: ", p); next }
    lab <- do.call(rbind, lapply(split(fit, fit$key), function(d) d[which.max(d$pt), ]))
    yr <- diff(range(df$expression, na.rm = TRUE))
    if(yr == 0) yr <- 1
    lab <- lab[order(lab$fit), ]
    ly <- lab$fit
    gap <- 0.055 * yr
    if(length(ly) > 1){
      for(i in 2:length(ly)) ly[i] <- max(ly[i], ly[i - 1] + gap)
    }
    lab$lab_y <- ly
    lab$x_at <- lab$pt
    ytop <- max(c(fit$fit, lab$lab_y), na.rm = TRUE) * 1.06
    res[[p]] <- ggplot() +
      geom_line(data = fit, aes(x = pt, y = fit, colour = gene), linewidth = 0.8) +
      geom_text(data = lab, aes(x = x_at, y = lab_y, label = key, colour = gene),
                hjust = 1.04, vjust = 0.5, size = label_size, fontface = "italic", show.legend = FALSE) +
      scale_colour_manual(values = gene_colors, guide = "none") +
      scale_x_continuous(expand = c(0.02, 0.02)) +
      scale_y_continuous(limits = c(0, ytop), expand = c(0, 0)) +
      labs(title = p, x = xlab, y = ylab) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9),
            axis.text = element_text(color = "black", size = 7))
  }
  res
}

resA <- plot_pt_gene_by_patient(pbmc_mono, pt_gene_hits, ncol = 3, ci = TRUE)
resA[["182LPIN2"]]
resA[["157LACC1"]]
resA[["183NOD2"]]

resB <- plot_pt_gene_overlay(pbmc_mono, pt_gene_hits, normalize = TRUE)
resB[["190PSTPIP1"]]
resB[["182LPIN2"]]
library(patchwork)
ord <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
resB_all <- wrap_plots(resB[ord], ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = NULL) &
  theme(plot.margin = margin(4, 4, 4, 4))
resB_all




