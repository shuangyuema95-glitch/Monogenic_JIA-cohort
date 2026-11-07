library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ggplotify)
library(cowplot)
library(patchwork)

sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6","157LACC1",
                  "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")

mono_metab_sub <- subset(pbmc1, celltype %in% c("CD14 Monocyte","CD16 Monocyte"))
ncol(mono_metab_sub);table(pbmc1$celltype)
Idents(mono_metab_sub) <- "celltype"
table(mono_metab_sub$celltype)
rm(pbmc1);gc()

countexp.Seurat <- sc.metabolism.Seurat(obj = mono_metab_sub,
                                        method = "VISION",
                                        imputation = F,
                                        ncores = 8,
                                        metabolism.type = "KEGG")

colnames(countexp.Seurat@assays$METABOLISM$score) <- gsub("\\.","-",colnames(countexp.Seurat@assays$METABOLISM$score))
score_change <- countexp.Seurat@assays$METABOLISM$score
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change))
metab_all_cols <- rownames(countexp.Seurat@assays$METABOLISM$score)
metab_score_df <- countexp.Seurat@meta.data[,c(metab_all_cols,"celltype","sample"),drop=FALSE]
unique(metab_score_df$sample)
colnames(metab_score_df)
#saveRDS(metab_score_df,file="./CellPreprocess/1105_scMetabolismScore.RDS")
metab_score_df<-readRDS("./CellPreprocess/1105_scMetabolismScore.RDS")


####(1)pheatmap for metabolism score
keep_pathways<-c("Glycolysis / Gluconeogenesis","Oxidative phosphorylation","Fatty acid biosynthesis","Fatty acid degradation","Glycerolipid metabolism","Glycerophospholipid metabolism","Ether lipid metabolism","Sphingolipid metabolism","Arachidonic acid metabolism","Pyruvate metabolism","Citrate cycle (TCA cycle)","Pentose phosphate pathway","Glutathione metabolism")
library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales)

sample_level <- c(paste0("C", 1:6), paste0("poly", 1:6),
                  "157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
mono_ids <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")

plot_pathway_heatmap_gg <- function(data, cells, pathways, sample_order = sample_level,
                                    pal = colorRampPalette(c("#542788","#F7F7F7","#D6604D"))(101),
                                    z_lim = 1.5, gap_width = 0.10,
                                    group_colors = c("HC" = "gray90","Polygenic" = "#F1BAAB","Monogenic" = "#D9C2D9"),
                                    cluster_rows = TRUE, clustering_method = "ward.D2",
                                    scale_mode = "row", border_col = NA, border_lwd = 0,
                                    row_labels = NULL, aspect = NULL,
                                    cell_width = 4, cell_height = 4,
                                    fontsize_row = 9, fontsize_col = 8, legend_title = "z-score"){
  gmap <- setNames(ifelse(sample_order %in% paste0("C", 1:6), "HC",
                          ifelse(sample_order %in% c(paste0("poly", 1:6), paste0("ploy", 1:6)), "Polygenic",
                                 ifelse(sample_order %in% mono_ids, "Monogenic", NA))), sample_order)
  out <- list()
  for(ct in cells){
    subdata <- data %>% filter(celltype == ct)
    avg <- subdata %>% pivot_longer(cols = all_of(pathways), names_to = "Pathway", values_to = "Score") %>%
      group_by(sample, Pathway) %>% summarise(v = mean(Score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = sample, values_from = v)
    mat <- avg %>% column_to_rownames("Pathway") %>% as.matrix()
    cols <- intersect(sample_order, colnames(mat))
    mat <- mat[, cols, drop = FALSE]
    if(scale_mode == "row") mat <- t(scale(t(mat)))
    row_ord <- rownames(mat)
    if(cluster_rows && nrow(mat) > 2){
      hc <- hclust(dist(mat), method = clustering_method)
      row_ord <- rownames(mat)[hc$order]
      mat <- mat[row_ord, , drop = FALSE]
    }
    y_lab <- if(!is.null(row_labels)){
      v <- row_labels[row_ord]; ifelse(is.na(v), row_ord, unname(v))
    } else row_ord
    grp <- unname(gmap[cols])
    xp <- numeric(length(cols)); shift <- 0
    for(k in seq_along(cols)){
      xp[k] <- k + shift
      if(k < length(cols) && grp[k+1] != grp[k]) shift <- shift + gap_width
    }
    names(xp) <- cols
    yp <- setNames(nrow(mat):1, row_ord)
    df <- data.frame(Pathway = rep(row_ord, times = ncol(mat)),
                     sample = rep(cols, each = nrow(mat)),
                     value = as.vector(mat), stringsAsFactors = FALSE)
    df$xpos <- as.numeric(xp[df$sample]); df$ypos <- as.numeric(yp[df$Pathway])
    n_row <- nrow(mat)
    asp <- if(!is.null(aspect)) aspect else ((n_row + 1) * cell_height) / (max(xp) * cell_width)
    p <- ggplot(df, aes(x = xpos, y = ypos, fill = value)) +
      geom_tile(width = 1, height = 1, color = border_col, linewidth = border_lwd) +
      scale_fill_gradientn(colors = pal, limits = c(-z_lim, z_lim), oob = squish, name = legend_title) +
      scale_x_continuous(breaks = as.numeric(xp), labels = cols, expand = c(0, 0)) +
      scale_y_continuous(breaks = n_row:1, labels = y_lab, expand = c(0, 0)) +
      coord_cartesian(xlim = c(0.5, max(xp) + 0.5), ylim = c(0.5, n_row + 1.5)) +
      labs(x = NULL, y = NULL, title = ct) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = fontsize_col),
            axis.text.y = element_text(color = "black", size = fontsize_row),
            axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9),
            legend.title = element_text(size = 8), legend.text = element_text(size = 8),
            legend.key.height = unit(4, "mm"), legend.key.width = unit(3, "mm"),
            aspect.ratio = asp,
            plot.margin = margin(4, 8, 4, 2))
    for(k in seq_along(cols)){
      p <- p + annotate("rect", xmin = xp[k] - 0.5, xmax = xp[k] + 0.5,
                        ymin = n_row + 0.5, ymax = n_row + 1.5,
                        fill = group_colors[grp[k]], color = NA)
    }
    attr(p, "fig_mm") <- c(width = max(xp) * cell_width + 42, height = (n_row + 1) * cell_height + 16)
    cat(sprintf("[%s] cell %.1f x %.1f mm -> fig %.0f x %.0f mm\n", ct, cell_width, cell_height,
                max(xp) * cell_width + 42, (n_row + 1) * cell_height + 16))
    out[[ct]] <- p
  }
  names(out) <- cells
  out
}

Metab_heats <- plot_pathway_heatmap_gg(data = metab_score_df,
                                       cells = c("CD14 Monocyte","CD16 Monocyte"),
                                       pathways = keep_pathways, row_labels = row_lab,
                                       z_lim = 1.5, gap_width = 0.10, border_col = NA,
                                       cell_width = 7.4, cell_height = 5.4)
wrap_plots(Metab_heats, ncol = 2)


####(2) error bar plot for metabolsim score
library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
Metab_Group_Bar <- function(data, pathways, cells,
                            celltype_levels = c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T","γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC","Platelet","Erythrocyte"),
                            hc_fill = "gray90", poly_fill = "#9CD2ED", mono_fill = "#D9C2D9",
                            hc_pt_col = "gray50", poly_pt_col = "#4E8FB8",
                            mono_cols = c("157LACC1" = "#8B9CD8","182LPIN2" = "#C47CA6","181NOD2" = "#A5BE76","183NOD2" = "#83BE62","200PSTPIP1" = "#5C1778","190PSTPIP1" = "#EFD97E"),
                            mono_shapes = c("157LACC1" = 16,"182LPIN2" = 17,"181NOD2" = 15,"183NOD2" = 18,"200PSTPIP1" = 8,"190PSTPIP1" = 7),
                            bar_width = 0.2, grp_off = 0.26, bar_gap = 0.32,
                            bar_lwd = 0.3, pt_size = 1.5, hc_poly_pt_size = 0.7,
                            pt_alpha = 0.95, jitter_w = 0.02, vline_col = "gray70", seed = 1,
                            y_pad = 0.12){
  mono_ids <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  DATA <- data[, c("celltype","sample", pathways)] %>%
    pivot_longer(cols = all_of(pathways), names_to = "pathway", values_to = "score")
  DATA$sample <- as.character(DATA$sample); DATA$celltype <- as.character(DATA$celltype)
  DATA <- DATA[!is.na(DATA$score) & DATA$celltype %in% cells, ]
  DATA$group <- ifelse(DATA$sample %in% paste0("C", 1:6), "HC",
                       ifelse(DATA$sample %in% c(paste0("poly", 1:6), paste0("ploy", 1:6)), "Polygenic",
                              ifelse(DATA$sample %in% mono_ids, "Monogenic", NA)))
  DATA <- DATA[!is.na(DATA$group), ]
  DATA$group <- factor(DATA$group, levels = c("HC","Polygenic","Monogenic"))
  DATA$pathway <- factor(DATA$pathway, levels = pathways)
  smry <- DATA %>% group_by(celltype, pathway, group, sample) %>% summarise(score = mean(score), .groups = "drop")
  smry$gene_label <- ifelse(smry$group == "Monogenic", smry$sample, NA)
  stat <- smry %>% group_by(celltype, pathway, group) %>%
    summarise(n = n(), mean_score = mean(score), median_score = median(score), .groups = "drop")
  step <- 1 + bar_gap
  pw_pos <- setNames((seq_along(pathways) - 1) * step + 1, pathways)
  centers <- as.numeric(pw_pos)
  bounds <- c(centers - step/2, centers[length(centers)] + step/2)
  off <- c("HC" = -grp_off, "Polygenic" = 0, "Monogenic" = grp_off)
  p_list <- list(); kw_list <- list(); pw_list <- list()
  for(ct in cells){
    d <- smry[smry$celltype == ct, ]; st <- stat[stat$celltype == ct, ]
    if(nrow(d) == 0) next
    d$xpos <- as.numeric(pw_pos[as.character(d$pathway)]) + as.numeric(off[as.character(d$group)])
    st$xpos <- as.numeric(pw_pos[as.character(st$pathway)]) + as.numeric(off[as.character(st$group)])
    y_hi <- max(st$mean_score, na.rm = TRUE) * (1 + y_pad)
    p <- ggplot() +
      geom_vline(xintercept = bounds, linetype = "dashed", color = vline_col, linewidth = 0.3) +
      geom_col(data = st, aes(x = xpos, y = mean_score, fill = group), width = bar_width, color = "black", linewidth = bar_lwd) +
      geom_point(data = d[d$group %in% c("HC","Polygenic"), , drop = FALSE],
                 aes(x = xpos, y = score, color = group), size = hc_poly_pt_size, alpha = pt_alpha,
                 position = position_jitter(width = jitter_w, height = 0, seed = seed), show.legend = FALSE) +
      geom_point(data = d[d$group == "Monogenic", , drop = FALSE],
                 aes(x = xpos, y = score, color = gene_label, shape = gene_label), size = pt_size, alpha = pt_alpha, stroke = 0.5,
                 position = position_jitter(width = jitter_w, height = 0, seed = seed), show.legend = FALSE) +
      scale_fill_manual(values = c("HC" = hc_fill, "Polygenic" = poly_fill, "Monogenic" = mono_fill)) +
      scale_color_manual(values = c("HC" = hc_pt_col, "Polygenic" = poly_pt_col, mono_cols)) +
      scale_shape_manual(values = mono_shapes) +
      scale_x_continuous(breaks = centers, labels = pathways, expand = expansion(add = 0.05)) +
      coord_cartesian(ylim = c(0, y_hi)) +
      annotate("text", x = 0.6, y = y_hi, hjust = 0, vjust = 1, size = 2.5, color = "black", label = ct) +
      labs(x = NULL, y = "Metabolism score") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
            axis.text.y = element_text(color = "black", size = 8),
            axis.title.y = element_text(color = "black", size = 9),
            legend.position = "none")
    p_list[[ct]] <- p
    pairs <- list(c("HC","Polygenic"), c("HC","Monogenic"), c("Polygenic","Monogenic"))
    for(pw in pathways){
      s <- smry[smry$celltype == ct & smry$pathway == pw, ]
      s$group <- droplevels(s$group)
      if(nlevels(s$group) < 2) next
      kw <- kruskal.test(score ~ group, data = s)
      kw_list[[length(kw_list)+1]] <- data.frame(celltype = ct, pathway = pw, n_cell = nrow(s),
                                                 n_group = nlevels(s$group), kw_p = kw$p.value, stringsAsFactors = FALSE)
      for(pp in pairs){
        v1 <- s$score[s$group == pp[1]]; v2 <- s$score[s$group == pp[2]]
        if(length(v1) < 2 || length(v2) < 2) next
        wt <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))
        pw_list[[length(pw_list)+1]] <- data.frame(celltype = ct, pathway = pw,
                                                   comparison = paste0(pp[1], " vs ", pp[2]), n1 = length(v1), n2 = length(v2),
                                                   mean1 = mean(v1), mean2 = mean(v2), median1 = median(v1), median2 = median(v2),
                                                   p_value = wt$p.value, stringsAsFactors = FALSE)
      }
    }
  }
  kw_tab <- do.call(rbind, kw_list); kw_tab$kw_padj <- p.adjust(kw_tab$kw_p, method = "BH")
  pair_tab <- do.call(rbind, pw_list); pair_tab$p_adj <- p.adjust(pair_tab$p_value, method = "BH")
  pair_tab$sig <- ifelse(is.na(pair_tab$p_adj), "ns", ifelse(pair_tab$p_adj < 0.001, "***",
                                                             ifelse(pair_tab$p_adj < 0.01, "**", ifelse(pair_tab$p_adj < 0.05, "*", "ns"))))
  list(plots = p_list,
       mean_tab = stat[, c("celltype","pathway","group","n","mean_score")],
       median_tab = stat[, c("celltype","pathway","group","n","median_score")],
       kw_tab = kw_tab, pair_tab = pair_tab)
}

bar_res <- Metab_Group_Bar(metab_score_df, pathways = keep_pathways, cells = c("CD14 Monocyte","CD16 Monocyte"),
                           bar_width = 0.2, grp_off = 0.26, bar_gap = 0.32,
                           bar_lwd = 0.3, pt_size = 1.5, hc_poly_pt_size = 0.7,
                           pt_alpha = 0.95, jitter_w = 0.02, seed = 1)
print(bar_res$plots[[1]])
print(bar_res$mean_tab, row.names = FALSE)
print(bar_res$median_tab, row.names = FALSE)
print(bar_res$kw_tab[bar_res$kw_tab$kw_padj < 0.05, ], row.names = FALSE)
print(bar_res$pair_tab[bar_res$pair_tab$p_adj < 0.05, ], row.names = FALSE)
wrap_plots(bar_res$plots, ncol = 1)


####(3) boxplot for scMetabolism score
library(dplyr); library(ggplot2)
Metab_Group_Box <- function(data, pathway, cells,
                            celltype_levels = c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T","γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC","Platelet","Erythrocyte"),
                            hc_col = "gray90", poly_col = "#9CD2ED", mono_col = "#D9C2D9",
                            hc_pt_col = "gray50", poly_pt_col = "#4E8FB8", box_col = "black",
                            mono_cols = c("157LACC1" = "#8B9CD8","182LPIN2" = "#C47CA6","181NOD2" = "#A5BE76","183NOD2" = "#83BE62","200PSTPIP1" = "#5C1778","190PSTPIP1" = "#EFD97E"),
                            mono_shapes = c("157LACC1" = 16,"182LPIN2" = 17,"181NOD2" = 15,"183NOD2" = 18,"200PSTPIP1" = 8,"190PSTPIP1" = 7),
                            cell_gap = 0.33, box_width = 0.27, grp_off = 0.28,
                            box_lwd = 0.3, box_fatten = 1.2, pt_size = 1.5,
                            hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02,
                            vline_col = "gray70", seed = 1){
  mono_ids <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  if(!"sample" %in% colnames(data)) data$sample <- data$datasets
  DATA <- data[, c(pathway, "celltype", "sample")]
  colnames(DATA) <- c("score","celltype","sample")
  DATA$sample <- as.character(DATA$sample); DATA$celltype <- as.character(DATA$celltype)
  DATA <- DATA[!is.na(DATA$score) & DATA$celltype %in% cells & !is.na(DATA$sample), ]
  DATA$celltype <- factor(DATA$celltype, levels = celltype_levels)
  DATA <- DATA[!is.na(DATA$celltype), ]
  DATA$group <- ifelse(DATA$sample %in% paste0("C", 1:6), "HC",
                       ifelse(DATA$sample %in% paste0("poly", 1:6), "Polygenic",
                              ifelse(DATA$sample %in% mono_ids, "Monogenic", NA)))
  DATA <- DATA[!is.na(DATA$group), ]
  DATA$group <- factor(DATA$group, levels = c("HC","Polygenic","Monogenic"))
  smry <- DATA %>% group_by(celltype, group, sample) %>% summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  smry$sample <- as.character(smry$sample)
  smry$gene_label <- ifelse(smry$group == "Monogenic", smry$sample, NA)
  step <- 1 + cell_gap
  ct_pos <- setNames(1 + (seq_along(celltype_levels) - 1) * step, celltype_levels)
  off <- c("HC" = -grp_off, "Polygenic" = 0, "Monogenic" = grp_off)
  smry$xpos <- as.numeric(ct_pos[as.character(smry$celltype)]) + as.numeric(off[as.character(smry$group)])
  cts_used <- celltype_levels[celltype_levels %in% smry$celltype]
  centers <- as.numeric(ct_pos[cts_used])
  bounds <- c(centers - step/2, centers[length(centers)] + step/2)
  p <- ggplot() +
    geom_vline(xintercept = bounds, linetype = "dashed", color = vline_col, linewidth = 0.3) +
    geom_boxplot(data = smry, aes(x = xpos, y = score, group = interaction(celltype, group), fill = group),
                 width = box_width, linewidth = box_lwd, fatten = box_fatten, outlier.shape = NA, color = box_col, show.legend = FALSE) +
    geom_point(data = smry[smry$group %in% c("HC","Polygenic"), , drop = FALSE],
               aes(x = xpos, y = score, color = group), size = hc_poly_pt_size, alpha = pt_alpha,
               position = position_jitter(width = jitter_w, height = 0, seed = seed), show.legend = FALSE) +
    geom_point(data = smry[smry$group == "Monogenic", , drop = FALSE],
               aes(x = xpos, y = score, color = gene_label, shape = gene_label), size = pt_size, alpha = pt_alpha, stroke = 0.5,
               position = position_jitter(width = jitter_w, height = 0, seed = seed), show.legend = FALSE) +
    scale_fill_manual(values = c("HC" = hc_col, "Polygenic" = poly_col, "Monogenic" = mono_col)) +
    scale_color_manual(values = c("HC" = hc_pt_col, "Polygenic" = poly_pt_col, mono_cols)) +
    scale_shape_manual(values = mono_shapes) +
    scale_x_continuous(breaks = centers, labels = cts_used, expand = expansion(add = 0.05)) +
    labs(x = NULL, y = pathway) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 9),
          legend.position = "none")
  pairs <- list(c("HC","Polygenic"), c("HC","Monogenic"), c("Polygenic","Monogenic"))
  rows <- list()
  for(ct in cts_used){
    sub <- smry[smry$celltype == ct, ]
    for(pp in pairs){
      v1 <- sub$score[sub$group == pp[1]]; v2 <- sub$score[sub$group == pp[2]]
      if(length(v1) < 2 || length(v2) < 2) next
      wt <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))
      rows[[length(rows)+1]] <- data.frame(pathway = pathway, celltype = ct, group1 = pp[1], group2 = pp[2],
                                           n1 = length(v1), n2 = length(v2), median1 = median(v1), median2 = median(v2),
                                           p_value = wt$p.value, stringsAsFactors = FALSE)
    }
  }
  sig <- do.call(rbind, rows)
  if(is.null(sig)) sig <- data.frame(pathway = pathway, celltype = character(), group1 = character(), group2 = character(),
                                     n1 = integer(), n2 = integer(), median1 = numeric(), median2 = numeric(), p_value = numeric(), p_adj = numeric(), sig = character())
  sig$p_adj <- p.adjust(sig$p_value, method = "BH")
  sig$sig <- ifelse(is.na(sig$p_adj), "ns", ifelse(sig$p_adj < 0.001, "***", ifelse(sig$p_adj < 0.01, "**", ifelse(sig$p_adj < 0.05, "*", "ns"))))
  list(plot = p, sig = sig)
}

keep_pathways <- c("Glycolysis / Gluconeogenesis","Oxidative phosphorylation","Fatty acid biosynthesis","Fatty acid degradation","Glycerolipid metabolism","Glycerophospholipid metabolism","Ether lipid metabolism","Sphingolipid metabolism","Arachidonic acid metabolism","Pyruvate metabolism","Citrate cycle (TCA cycle)","Pentose phosphate pathway","Glutathione metabolism")
stopifnot(all(keep_pathways %in% colnames(metab_score_df)))

metab_res <- lapply(keep_pathways, function(pw) Metab_Group_Box(metab_score_df, pathway = pw, cells = c("CD14 Monocyte","CD16 Monocyte"), cell_gap = 0.33,
                                                                box_width = 0.27, grp_off = 0.28, box_lwd = 0.3,
                                                                box_fatten = 1.2, pt_size = 1.5,
                                                                hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02, seed = 1))
names(metab_res) <- keep_pathways
metab_sig <- do.call(rbind, lapply(keep_pathways, function(x) metab_res[[x]]$sig))
metab_sig$p_adj_global <- p.adjust(metab_sig$p_value, method = "BH")
print(metab_sig[metab_sig$sig != "ns", ], row.names = FALSE)

print(metab_res[["Glycolysis / Gluconeogenesis"]]$plot)
print(metab_res[["Arachidonic acid metabolism"]]$plot)
print(metab_res[["Oxidative phosphorylation"]]$plot)

plots_metab <- lapply(metab_res, function(x) x$plot)
wrap_plots(plots_metab, ncol = 4)
wrap_plots(lapply(plots_metab, function(p) p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())), nrow = 1)

