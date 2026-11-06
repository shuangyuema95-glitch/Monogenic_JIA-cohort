library(Seurat); library(data.table); library(ggplot2); library(dplyr)

######1. cytotrace2 score loading
cyto2 <- as.data.frame(fread("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\cytotrace2_results.txt"))
pbmc1 <- readRDS("E:\\Cohort PPT\\JIA\\code\\pbmc1.RDS")

pbmc1@meta.data$cytoscore <- cyto2[match(rownames(pbmc1@meta.data), cyto2$V1), "CytoTRACE2_Score"]
cat("Matched:", sum(!is.na(pbmc1@meta.data$cytoscore)), "/", ncol(pbmc1), "\n")
cat("Score range:", round(range(pbmc1@meta.data$cytoscore, na.rm = TRUE), 4), "\n")
group_levels <- c("HC","Polygenic","157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
mono_ids <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
group_colors <- c("HC"="#6699CC","Polygenic"="#FEC260","157LACC1"="#B5BBE3","182LPIN2"="#E2A9C9",
                  "181NOD2"="#CBDAA9","183NOD2"="#B1DA99","200PSTPIP1"="#7B1FA2","190PSTPIP1"="#FFEFC1")
celltype_levels <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T","γδT",
                     "Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC","Platelet","Erythrocyte")

pbmc1@meta.data$group <- NA
pbmc1@meta.data$group[pbmc1@meta.data$datasets %in% paste0("C", 1:6)] <- "HC"
pbmc1@meta.data$group[pbmc1@meta.data$datasets %in% c(paste0("poly", 1:6), paste0("ploy", 1:6))] <- "Polygenic"
pbmc1@meta.data$group[pbmc1@meta.data$datasets %in% mono_ids] <- pbmc1@meta.data$datasets[pbmc1@meta.data$datasets %in% mono_ids]
pbmc1@meta.data$group <- factor(pbmc1@meta.data$group, levels = group_levels)

meta_df <- pbmc1@meta.data[, c("celltype","group","cytoscore")]
meta_df <- meta_df[!is.na(meta_df$cytoscore) & !is.na(meta_df$group), ]
meta_df$celltype <- factor(meta_df$celltype, levels = celltype_levels)
meta_df <- meta_df[!is.na(meta_df$celltype), ]
all_ct <- levels(droplevels(meta_df$celltype))
cat("\n=== Celltypes:", length(all_ct), "===\n"); print(all_ct)

##### testing: Kruskal-Wallis + pairwise Wilcoxon (Mann-Whitney)
kw_list <- list(); pw_list <- list()
for(ct in all_ct){
  sub <- meta_df[meta_df$celltype == ct, ]
  sub$group <- droplevels(sub$group)
  grp_n <- table(sub$group)
  valid_grps <- names(grp_n[grp_n >= 30])
  sub <- sub[sub$group %in% valid_grps, ]
  sub$group <- droplevels(sub$group)
  if(nlevels(sub$group) < 2 || nrow(sub) < 10) next
  kw <- kruskal.test(cytoscore ~ group, data = sub)
  kw_list[[ct]] <- data.frame(celltype = ct, n_cell = nrow(sub), n_group = nlevels(sub$group),
                              p_kw = kw$p.value, stringsAsFactors = FALSE)
  pw <- pairwise.wilcox.test(sub$cytoscore, sub$group, p.adjust.method = "BH")
  pw_tab <- as.data.frame(as.table(pw$p.value))
  pw_tab <- pw_tab[!is.na(pw_tab$Freq), ]
  pw_list[[ct]] <- data.frame(celltype = ct, comparison = paste0(pw_tab$Var1, " vs ", pw_tab$Var2),
                              p_within = pw_tab$Freq, stringsAsFactors = FALSE)
}
kw_df <- do.call(rbind, kw_list); kw_df$p_global <- p.adjust(kw_df$p_kw, method = "BH")
pw_df <- do.call(rbind, pw_list); pw_df$p_global <- p.adjust(pw_df$p_within, method = "BH")
print(kw_df[order(kw_df$p_global), ], row.names = FALSE)
print(pw_df[pw_df$p_global < 0.05, ], row.names = FALSE)
# openxlsx::write.xlsx(list(KW = kw_df, PW = pw_df), "E:/Cohort PPT/JIA/code/CellTrajectory/Cytotrace2_test.xlsx")

##### plot
plot_cytotrace2 <- function(meta_df, kw_df, group_colors, celltype_levels, mode = "color",
                            box_width = 0.8, box_lwd = 0.8,
                            hlines = c("HC", "Polygenic"), hline_stat = "median", hline_lwd = 0.6,
                            x_labels = c("Polygenic" = "Polygenic JIA", "157LACC1" = "P7", "182LPIN2" = "P10",
                                         "181NOD2" = "P14", "183NOD2" = "P15", "200PSTPIP1" = "P17", "190PSTPIP1" = "P18"),
                            y_pad = 0.1, y_title = "CytoTRACE2 Score",
                            pt_col = "cytoscore", group_col = "group", celltype_col = "celltype",
                            group_order = c("HC","Polygenic","157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")){
  df <- meta_df[, c(celltype_col, group_col, pt_col)]
  colnames(df) <- c("celltype","group","score")
  df <- df[!is.na(df$score) & !is.na(df$group), ]
  df$group <- factor(as.character(df$group), levels = group_order)
  df <- df[!is.na(df$group), ]
  df$group_idx <- as.numeric(df$group)
  group_levels <- levels(df$group)
  display_labels <- group_levels
  if(!is.null(x_labels)) for(nm in names(x_labels)) if(nm %in% display_labels) display_labels[display_labels == nm] <- x_labels[nm]
  all_ct2 <- celltype_levels[celltype_levels %in% unique(as.character(df$celltype))]
  p_list <- list(); sig_celltypes <- c()
  for(ct in all_ct2){
    sub <- df[as.character(df$celltype) == ct, ]
    kw_row <- kw_df[kw_df$celltype == ct, ]
    if(nrow(kw_row) > 0){
      p_adj <- kw_row$p_global[1]
      lab_p <- ifelse(is.na(p_adj), "n/a", formatC(p_adj, format = "e", digits = 2))
      if(!is.na(p_adj) && p_adj < 0.05) sig_celltypes <- c(sig_celltypes, ct)
    } else lab_p <- "n/a (too few cells)"
    whisker_values <- c()
    for(g in levels(sub$group)){
      tmp <- sub$score[sub$group == g]
      if(length(tmp) > 1){ bp <- boxplot.stats(tmp); whisker_values <- c(whisker_values, bp$stats[1], bp$stats[5]) }
    }
    y_min <- min(whisker_values, na.rm = TRUE); y_max <- max(whisker_values, na.rm = TRUE)
    y_range <- y_max - y_min; if(y_range == 0) y_range <- 0.01
    y_lo <- y_min - y_range * 0.02; y_hi <- y_max + y_range * y_pad
    if(mode == "fill"){
      p <- ggplot(sub, aes(x = group_idx, y = score, fill = group)) +
        geom_boxplot(width = box_width, linewidth = box_lwd, color = "black", outlier.shape = NA, fatten = 1.2) +
        scale_fill_manual(values = group_colors)
    } else {
      p <- ggplot(sub, aes(x = group_idx, y = score)) +
        geom_boxplot(aes(color = group), width = box_width, linewidth = box_lwd, outlier.shape = NA, fill = "white", fatten = 1.2) +
        scale_color_manual(values = group_colors)
    }
    for(hg in hlines){
      if(hg %in% group_levels){
        hv <- if(hline_stat == "median") median(sub$score[sub$group == hg], na.rm = TRUE) else mean(sub$score[sub$group == hg], na.rm = TRUE)
        p <- p + geom_hline(yintercept = hv, linetype = "dashed", color = group_colors[hg], linewidth = hline_lwd, alpha = 0.9)
      }
    }
    p <- p + scale_x_continuous(breaks = 1:length(group_levels), labels = display_labels) +
      coord_cartesian(ylim = c(y_lo, y_hi)) +
      annotate("text", x = 0.66, y = y_hi, hjust = 0, vjust = 1, size = 2.5, color = "black",
               label = paste0(ct, "\nKW p.adj = ", lab_p)) +
      labs(x = NULL, y = y_title) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.3),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
            axis.text.y = element_text(color = "black", size = 9),
            axis.title = element_text(color = "black", size = 10),
            plot.margin = margin(2, 2, 2, 2), legend.position = "none")
    p_list[[ct]] <- p
  }
  cat("Mode:", mode, "| Cells:", length(p_list), "| KW sig:", paste(sig_celltypes, collapse = ", "), "\n")
  list(plots = p_list, significant = sig_celltypes)
}

res_color <- plot_cytotrace2(meta_df, kw_df, group_colors, celltype_levels = celltype_levels, mode = "color",
                             box_width = 0.8, box_lwd = 0.8, hline_stat = "median", hline_lwd = 0.6)
res_fill  <- plot_cytotrace2(meta_df, kw_df, group_colors, celltype_levels = celltype_levels, mode = "fill",
                             box_width = 0.8, box_lwd = 0.8, hline_stat = "median", hline_lwd = 0.6)

wrap_plots(lapply(res_color$plots, function(p) p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())), nrow = 1)
wrap_plots(res_color$plots, ncol = 1)
wrap_plots(lapply(res_fill$plots, function(p) p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())), nrow = 1)
wrap_plots(res_fill$plots, ncol = 1)

print(res_color$plots$`CD14 Monocyte`)
print(res_fill$plots$`CD14 Monocyte`)