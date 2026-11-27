library(Seurat); library(data.table); library(ggplot2); library(dplyr); library(patchwork)

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
mono_ct <- c("CD14 Monocyte","CD16 Monocyte")
x_labels_def <- c("Polygenic"="Polygenic JIA","157LACC1"="P7","182LPIN2"="P10",
                  "181NOD2"="P14","183NOD2"="P15","200PSTPIP1"="P17","190PSTPIP1"="P18")

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

######2. testing per celltype: Kruskal-Wallis + pairwise Wilcoxon
kw_list <- list(); pw_list <- list()
for(ct in all_ct){
  sub <- meta_df[meta_df$celltype == ct, ]
  sub$group <- droplevels(sub$group)
  grp_n <- table(sub$group)
  sub <- sub[sub$group %in% names(grp_n[grp_n >= 30]), ]
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

######3. testing pooled monocyte (CD14 + CD16 merged)
mono_df <- meta_df[as.character(meta_df$celltype) %in% mono_ct, c("celltype","group","cytoscore")]
mono_df$celltype <- "Monocyte"
mono_df$group <- droplevels(mono_df$group)
grp_n <- table(mono_df$group)
mono_df <- mono_df[mono_df$group %in% names(grp_n[grp_n >= 30]), ]
mono_df$group <- droplevels(mono_df$group)
cat("\n=== Monocyte cells:", nrow(mono_df), "| groups:", nlevels(mono_df$group), "===\n")
print(table(mono_df$group))

kw_mono <- kruskal.test(cytoscore ~ group, data = mono_df)
kw_mono_df <- data.frame(celltype = "Monocyte", n_cell = nrow(mono_df),
                         n_group = nlevels(mono_df$group), p_kw = kw_mono$p.value,
                         p_global = NA_real_, stringsAsFactors = FALSE)
pw_mono <- pairwise.wilcox.test(mono_df$cytoscore, mono_df$group, p.adjust.method = "BH")
pw_tab <- as.data.frame(as.table(pw_mono$p.value))
pw_tab <- pw_tab[!is.na(pw_tab$Freq), ]
pw_mono_df <- data.frame(celltype = "Monocyte", comparison = paste0(pw_tab$Var1, " vs ", pw_tab$Var2),
                         p_within = pw_tab$Freq, p_global = NA_real_, stringsAsFactors = FALSE)

kw_all <- rbind(kw_df, kw_mono_df); kw_all$p_global <- p.adjust(kw_all$p_kw, method = "BH")
pw_all <- rbind(pw_df, pw_mono_df); pw_all$p_global <- p.adjust(pw_all$p_within, method = "BH")
kw_mono_df$p_global <- kw_all$p_global[kw_all$celltype == "Monocyte"]
pw_mono_df$p_global <- pw_all$p_global[pw_all$celltype == "Monocyte"]
print(kw_all[kw_all$celltype == "Monocyte", ], row.names = FALSE)
print(pw_all[pw_all$celltype == "Monocyte" & pw_all$p_global < 0.05, ], row.names = FALSE)

######4. plot per celltype (median / mean switchable)
plot_cytotrace2 <- function(meta_df, kw_df, group_colors, celltype_levels, mode = "color",
                            box_width = 0.8, box_lwd = 0.8,
                            hlines = c("HC","Polygenic"), hline_stat = "median", hline_lwd = 0.6,
                            hline_lab = TRUE, hline_label_size = 2.2,
                            x_labels = x_labels_def, y_pad = 0.1, y_title = "CytoTRACE2 Score",
                            pt_col = "cytoscore", group_col = "group", celltype_col = "celltype",
                            group_order = group_levels){
  hline_stat <- match.arg(hline_stat, c("median","mean"))
  mode <- match.arg(mode, c("color","fill"))
  df <- meta_df[, c(celltype_col, group_col, pt_col)]
  colnames(df) <- c("celltype","group","score")
  df <- df[!is.na(df$score) & !is.na(df$group), ]
  df$group <- factor(as.character(df$group), levels = group_order)
  df <- df[!is.na(df$group), ]
  df$group_idx <- as.numeric(df$group)
  group_levels2 <- levels(df$group)
  display_labels <- group_levels2
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
      if(hg %in% group_levels2){
        hv <- if(hline_stat == "median") median(sub$score[sub$group == hg], na.rm = TRUE) else mean(sub$score[sub$group == hg], na.rm = TRUE)
        p <- p + geom_hline(yintercept = hv, linetype = "dashed", color = group_colors[hg], linewidth = hline_lwd, alpha = 0.9)
        if(hline_lab){
          hl_x <- if(hg == "HC") 0.72 else length(group_levels2) + 0.35
          p <- p + annotate("text", x = hl_x, y = hv, hjust = 0, vjust = -0.4, size = hline_label_size,
                            color = group_colors[hg],
                            label = paste0(hg, " ", hline_stat, " = ",
                                           formatC(hv, format = "f", digits = 3, drop0trailing = TRUE)))
        }
      }
    }
    p <- p + scale_x_continuous(breaks = 1:length(group_levels2), labels = display_labels) +
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
  cat("Mode:", mode, "| hline:", hline_stat, "| panels:", length(p_list), "| KW sig:",
      paste(sig_celltypes, collapse = ", "), "\n")
  list(plots = p_list, significant = sig_celltypes)
}

######5. plot pooled monocyte (CD14 + CD16 merged), median / mean switchable
plot_mono_pooled <- function(meta_df, kw_mono_df = NULL, group_colors, group_order,
                             cts = mono_ct, box_width = 0.8, box_lwd = 0.8,
                             hlines = c("HC","Polygenic"), hline_stat = "median", hline_lwd = 0.6,
                             hline_lab = TRUE, hline_label_size = 2.2,
                             x_labels = x_labels_def, y_pad = 0.1,
                             y_title = "CytoTRACE2 Score", title = "Monocyte"){
  hline_stat <- match.arg(hline_stat, c("median","mean"))
  df <- meta_df[as.character(meta_df$celltype) %in% cts, c("celltype","group","cytoscore")]
  colnames(df) <- c("celltype","group","score")
  df <- df[!is.na(df$score) & !is.na(df$group), ]
  df$celltype <- "Monocyte"
  gl <- group_order[group_order %in% as.character(df$group)]
  df$group <- factor(as.character(df$group), levels = gl)
  df <- df[!is.na(df$group), ]
  df$group_idx <- as.numeric(df$group)
  
  if(is.null(kw_mono_df)){
    sub <- df; sub$group <- droplevels(sub$group)
    p_adj <- if(nlevels(sub$group) >= 2 && nrow(sub) >= 10) kruskal.test(score ~ group, data = sub)$p.value else NA
  } else {
    p_adj <- kw_mono_df$p_global[1]
  }
  lab_p <- ifelse(is.na(p_adj), "n/a", formatC(p_adj, format = "e", digits = 2))
  
  lab <- gl
  for(nm in names(x_labels)) if(nm %in% lab) lab[lab == nm] <- x_labels[nm]
  
  w <- c()
  for(g in gl){
    t <- df$score[df$group == g]
    if(length(t) > 1){ b <- boxplot.stats(t); w <- c(w, b$stats[1], b$stats[5]) }
  }
  if(length(w) == 0) w <- range(df$score, na.rm = TRUE)
  yr <- diff(range(w)); if(yr == 0) yr <- 0.01
  y_lo <- min(w) - yr * 0.02; y_hi <- max(w) + yr * y_pad
  
  p <- ggplot(df, aes(x = group_idx, y = score, fill = group)) +
    geom_boxplot(width = box_width, linewidth = box_lwd, color = "black", outlier.shape = NA, fatten = 1.2) +
    scale_fill_manual(values = group_colors)
  for(hg in hlines){
    if(hg %in% gl){
      hv <- if(hline_stat == "median") median(df$score[df$group == hg], na.rm = TRUE) else mean(df$score[df$group == hg], na.rm = TRUE)
      p <- p + geom_hline(yintercept = hv, linetype = "dashed", color = group_colors[hg], linewidth = hline_lwd, alpha = 0.9)
      if(hline_lab){
        hl_x <- if(hg == "HC") 0.72 else length(gl) + 0.35
        p <- p + annotate("text", x = hl_x, y = hv, hjust = 0, vjust = -0.4, size = hline_label_size,
                          color = group_colors[hg],
                          label = paste0(hg, " ", hline_stat, " = ",
                                         formatC(hv, format = "f", digits = 3, drop0trailing = TRUE)))
      }
    }
  }
  p + scale_x_continuous(breaks = seq_along(gl), labels = lab) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    annotate("text", x = 0.66, y = y_hi, hjust = 0, vjust = 1, size = 2.5, color = "black",
             label = paste0(title, "\nKW p.adj = ", lab_p)) +
    labs(x = NULL, y = y_title, title = title) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          axis.line = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.3),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
          axis.text.y = element_text(color = "black", size = 9),
          axis.title = element_text(color = "black", size = 10),
          plot.title = element_text(color = "black", size = 11, hjust = 0),
          plot.margin = margin(2, 2, 2, 2), legend.position = "none")
}

######6. run
res_color_med <- plot_cytotrace2(meta_df, kw_df, group_colors, celltype_levels = celltype_levels,
                                 mode = "color", hline_stat = "median")
res_fill_med  <- plot_cytotrace2(meta_df, kw_df, group_colors, celltype_levels = celltype_levels,
                                 mode = "fill", hline_stat = "median")
res_fill_mean <- plot_cytotrace2(meta_df, kw_df, group_colors, celltype_levels = celltype_levels,
                                 mode = "fill", hline_stat = "mean")

mono_med  <- plot_mono_pooled(meta_df, kw_mono_df, group_colors, group_order = group_levels, hline_stat = "median")
mono_mean <- plot_mono_pooled(meta_df, kw_mono_df, group_colors, group_order = group_levels, hline_stat = "mean")

wrap_plots(res_fill_med$plots, nrow = 1)
wrap_plots(res_fill_med$plots, ncol = 1)
p1 <- res_fill_med$plots$`CD14 Monocyte`; p2 <- res_fill_med$plots$`CD16 Monocyte`; p3 <- res_fill_med$plots$pDC
p1 / p2 / p3
print(mono_med); print(mono_mean)

######7. export test results
library(openxlsx)
out_dir <- "E:/Cohort PPT/JIA/code/CellTrajectory"
out_file <- file.path(out_dir, paste0("CytoTRACE2_KW_Wilcoxon_celltype_and_pooledMonocyte_", Sys.Date(), ".xlsx"))
##KW kruskail-wails; PW: MU
sheet_list <- list(
  KW_celltype            = kw_df[order(kw_df$p_global), ],
  PW_celltype            = pw_df[order(pw_df$p_global), ],
  KW_pooledMonocyte      = kw_mono_df,
  PW_pooledMonocyte      = pw_mono_df[order(pw_mono_df$p_global), ],
  KW_all_with_monocyte   = kw_all[order(kw_all$p_global), ],
  PW_all_with_monocyte   = pw_all[order(pw_all$p_global), ]
)

write.xlsx(sheet_list, file = out_file, overwrite = TRUE)
cat("Saved:", out_file, "\n")

