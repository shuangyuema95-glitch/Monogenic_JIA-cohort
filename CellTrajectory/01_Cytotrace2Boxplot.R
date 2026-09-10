
library(Seurat); library(data.table); library(ggplot2); library(dplyr)

######1. cytotrace2 score loading
cyto2 <- as.data.frame(fread("E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\cytotrace2_results.txt"))
pbmc1 <- readRDS("E:\\Cohort PPT\\JIA\\code\\pbmc1.RDS")

pbmc1@meta.data$cytoscore <- cyto2[match(rownames(pbmc1@meta.data), cyto2$V1), "CytoTRACE2_Score"]
cat("Matched:", sum(!is.na(pbmc1@meta.data$cytoscore)), "/", ncol(pbmc1), "\n")
cat("Score range:", round(range(pbmc1@meta.data$cytoscore, na.rm = TRUE), 4), "\n")

group_levels <- c("HC","Polygenic","181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
pbmc1@meta.data$group <- NA
pbmc1@meta.data$group[pbmc1@meta.data$datasets %in% paste0("C", 1:6)] <- "HC"
pbmc1@meta.data$group[pbmc1@meta.data$datasets %in% paste0("ploy", 1:6)] <- "Polygenic"
pbmc1@meta.data$group[pbmc1@meta.data$datasets == "181NOD2"]   <- "181NOD2"
pbmc1@meta.data$group[pbmc1@meta.data$datasets == "182LPIN2"]   <- "182LPIN2"
pbmc1@meta.data$group[pbmc1@meta.data$datasets == "183NOD2"]   <- "183NOD2"
pbmc1@meta.data$group[pbmc1@meta.data$datasets == "190PSTPIP1"] <- "190PSTPIP1"
pbmc1@meta.data$group[pbmc1@meta.data$datasets == "200PSTPIP1"] <- "200PSTPIP1"
pbmc1@meta.data$group <- factor(pbmc1@meta.data$group, levels = group_levels)


meta_df <- pbmc1@meta.data[, c("celltype","group","cytoscore")]
meta_df <- meta_df[!is.na(meta_df$cytoscore) & !is.na(meta_df$group), ]
all_ct <- sort(unique(meta_df$celltype))
cat("\n=== Celltypes:", length(all_ct), "===\n")
print(all_ct)

######2. testing
test_list <- list()
for (ct in all_ct) {
  sub <- meta_df[meta_df$celltype == ct, ]
  grp_n <- table(sub$group)
  valid_grps <- names(grp_n[grp_n >= 30])
  sub <- sub[sub$group %in% valid_grps, ]
  if (length(unique(sub$group)) < 2) next
  
  kw <- kruskal.test(cytoscore ~ group, data = sub)
  test_list[[length(test_list) + 1]] <- data.frame(
    celltype = ct, comparison = "Overall (Kruskal-Wallis)",
    p_value = kw$p.value, stringsAsFactors = FALSE)
  
  pw <- pairwise.wilcox.test(sub$cytoscore, sub$group, p.adjust.method = "BH")
  pw_df <- as.data.frame(as.table(pw$p.value))
  pw_df <- pw_df[!is.na(pw_df$Freq), ]
  for (i in 1:nrow(pw_df)) {
    test_list[[length(test_list) + 1]] <- data.frame(
      celltype = ct,
      comparison = paste0(pw_df$Var1[i], " vs ", pw_df$Var2[i]),
      p_value = pw_df$Freq[i], stringsAsFactors = FALSE)
  }
}
test_df <- do.call(rbind, test_list)
test_df$p_adj <- p.adjust(test_df$p_value, method = "BH")
test_df <- test_df[order(test_df$celltype, test_df$p_adj), ]

cat("\n=== Significant overall tests (p_adj < 0.05) ===\n")
print(test_df[test_df$comparison == "Overall (Kruskal-Wallis)" & test_df$p_adj < 0.05, ])
openxlsx::write.xlsx(test_df,file="E:\\Cohort PPT\\JIA\\code\\CellTrajectory\\Cytotrace2_test.xlsx")

######3 plot as boxplots
plot_cytotrace2 <- function(meta_df, test_df, group_colors, mode = "color",
                            box_width = 0.6, box_lwd = 0.5,
                            loess_color = "black", loess_lwd = 0.8,
                            x_labels = c("182LPIN2"="P9","181NOD2"="P13","183NOD2"="P14",
                                         "200PSTPIP1"="P15","190PSTPIP1"="P16"),
                            y_pad = 0.05, pt_col = "cytoscore",
                            group_col = "group", celltype_col = "celltype") {
  library(ggplot2); library(dplyr)
  df <- meta_df[, c(celltype_col, group_col, pt_col)]
  colnames(df) <- c("celltype","group","score")
  df <- df[!is.na(df$score) & !is.na(df$group), ]
  group_order <- c("HC","Polygenic","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  df$group <- factor(df$group, levels = group_order)
  df <- df[!is.na(df$group), ]
  df$group_idx <- as.numeric(df$group)
  group_levels <- levels(df$group)
  display_labels <- group_levels
  if (!is.null(x_labels)) {
    for (nm in names(x_labels)) if (nm %in% display_labels) display_labels[display_labels == nm] <- x_labels[nm]
  }
  all_ct <- sort(unique(df$celltype))
  p_list <- list(); sig_celltypes <- c()
  for (ct in all_ct) {
    sub <- df[df$celltype == ct, ]
    kw_row <- test_df[test_df$celltype == ct & test_df$comparison == "Overall (Kruskal-Wallis)", ]
    if (nrow(kw_row) > 0) {
      p_adj <- kw_row$p_adj[1]
      p_label <- paste0("Kruskal-Wallis test\np.adj = ", formatC(p_adj, format = "e", digits = 2))
      if (p_adj < 0.05) sig_celltypes <- c(sig_celltypes, ct)
    } else { p_label <- "n/a (too few cells)" }
    whisker_values <- c()
    for (g in levels(sub$group)) {
      tmp <- sub$score[sub$group == g]
      if (length(tmp) > 1) { bp <- boxplot.stats(tmp); whisker_values <- c(whisker_values, bp$stats[1], bp$stats[5]) }
    }
    y_min <- min(whisker_values, na.rm = TRUE); y_max <- max(whisker_values, na.rm = TRUE)
    y_range <- y_max - y_min; if (y_range == 0) y_range <- 0.01
    y_lo <- y_min - y_range * y_pad; y_hi <- y_max + y_range * y_pad
    p <- ggplot(sub, aes(x = group_idx, y = score)) +
      geom_boxplot(aes(color = group), width = box_width, linewidth = box_lwd, outlier.shape = NA, fill = "white") +
      scale_color_manual(values = group_colors) +
      scale_x_continuous(breaks = 1:length(group_levels), labels = display_labels) +
      coord_cartesian(ylim = c(y_lo, y_hi)) +
      labs(x = NULL, y = "CytoTRACE2 Score", title = paste0(ct, "\n", p_label)) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.3),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
            axis.text.y = element_text(color = "black", size = 9),
            axis.title = element_text(color = "black", size = 10),
            plot.title = element_text(hjust = 0.5, size = 9), legend.position = "none")
    if (mode == "color2") {
      p <- p + geom_smooth(method = "loess", se = FALSE, color = loess_color, linewidth = loess_lwd)
    } else if (mode == "color3") {
      if ("HC" %in% group_levels) {
        hc_m <- mean(sub$score[sub$group == "HC"], na.rm = TRUE)
        p <- p + geom_hline(yintercept = hc_m, linetype = "dashed", color = group_colors["HC"], linewidth = 0.6)
      }
      if ("Polygenic" %in% group_levels) {
        poly_m <- mean(sub$score[sub$group == "Polygenic"], na.rm = TRUE)
        p <- p + geom_hline(yintercept = poly_m, linetype = "dashed", color = group_colors["Polygenic"], linewidth = 0.6)
      }
    }
    p_list[[ct]] <- p
  }
  cat("Mode:", mode, "| Significant:", paste(sig_celltypes, collapse = ", "), "\n")
  return(list(plots = p_list, significant = sig_celltypes))
}


# res_color <- plot_cytotrace2(meta_df, test_df, group_colors, mode = "color",
#                              box_width = 0.7, box_lwd = 0.8)
# res_c2    <- plot_cytotrace2(meta_df, test_df, group_colors, mode = "color2",
#                              box_lwd = 0.8, loess_color = "#333333", loess_lwd = 1)
res_c3<-plot_cytotrace2(meta_df, test_df, group_colors, mode = "color3",box_width = 0.8, box_lwd = 0.8)



library(patchwork)
library(patchwork)

p1 <- res_c3$plots$`CD14 Monocyte` +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

p2 <- res_c3$plots$`CD16 Monocyte` +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank())

p3 <- res_c3$plots$Plasma +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank())

p1 | p2 | p3
