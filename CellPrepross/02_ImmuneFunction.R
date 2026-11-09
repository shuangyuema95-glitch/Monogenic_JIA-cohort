library(Seurat)
library(patchwork)
library(ggplot2)
library(AUCell)
library(ggsignif)
library(ggpubr)
library(ggplotify)
library(msigdbr)
library(dplyr)

setwd("E:/Cohort PPT/JIA/code")
pbmc1<-readRDS("pbmc1.RDS")

#####1 ----Immune signature score by AUCell algorithm-----
###(1)score calculation
DefaultAssay(pbmc1) <- 'RNA'
score_Calu <- function(pbmc1, inter_gene, type) {
  library(AUCell)
  name <- paste0("auc_", type)
  UC <- AUCell_run(GetAssayData(pbmc1, assay = "RNA", slot = "counts"), inter_gene)
  identical(colnames(getAUC(UC)), rownames(pbmc1@meta.data))
  pbmc1$aucScore <- as.numeric(getAUC(UC))
  colnames(pbmc1@meta.data)[match("aucScore", colnames(pbmc1@meta.data))] <- name
  
  return(pbmc1)
}
load("E:\\通路基因集合\\gene4pathway.Rdata")
IFN28 <- c("DDX60","EPSTI1","HERC5","HERC6","IFI27","IFI44","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFIT5","ISG15","LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SIGLEC1","SPATS2L","USP18","CXCL10","GBP1","SOCS1")
pyrop <- c("AIM2","APIP","CASP1","CASP4","CASP8","DHX9","ELANE","GSDMA","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB","NAIP","NLRC4","NLRP1","NLRP6","NLRP9","TREM2","ZBP1","NLRP3","PYCARD","NEK7","IL1B")
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
auto_genes <- reactome %>% filter(gs_name == "REACTOME_AUTOPHAGY") %>% pull(gene_symbol) %>% unique()
hm <- msigdbr(species = "Homo sapiens", category = "H")
gly_genes <- hm %>% filter(gs_name == "HALLMARK_GLYCOLYSIS") %>% pull(gene_symbol) %>% unique()
oxphos_genes <- hm %>% filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% pull(gene_symbol) %>% unique()
arachidonic_genes <- c("PLA2G4A","PLA2G6","PLA2G2A","PTGS1","PTGS2","ALOX5","ALOX5AP","ALOX12","ALOX15","ALOX15B","LTC4S","PTGES","PTGES2","TBXAS1","GPX4","CYP2J2","CYP4F2","EPHX2","HPGDS")
inflam_genes <- hm %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique()


pbmc1 <- score_Calu(pbmc1, NFKB$gene, "NFKB")
pbmc1 <- score_Calu(pbmc1, IFN28, "IFNI")
pbmc1 <- score_Calu(pbmc1, IFN2, "IFNII")
pbmc1 <- score_Calu(pbmc1, MAPK$gene, "MAPK")
pbmc1 <- score_Calu(pbmc1, pyrop, "pyrop")
pbmc1 <- score_Calu(pbmc1, auto_genes, "autophagy")
pbmc1 <- score_Calu(pbmc1, gly_genes, "glycolysis")
pbmc1 <- score_Calu(pbmc1, oxphos_genes, "OXPHOS")
pbmc1 <- score_Calu(pbmc1, arachidonic_genes, "ARA")
pbmc1 <- score_Calu(pbmc1, inflam_genes, "Inflammatory")
scores<-pbmc1@meta.data
saveRDS(scores,file="1102_pbmc1_aucScore.rds")

###(2)boxplot visualization, by three groups
Violin_Score_Adaptive_Sig <- function(pbmc1, plot_mode = c("cell_facet","score_single")){
  plot_mode <- match.arg(plot_mode)
  DATA <- pbmc1@meta.data
  status_level <- c("Healthy controls","Polygenic JIA","Monogenic JIA")
  cell_color_vec3 <- c("#6699CC","#FEC260","#AF478A")
  color_map <- setNames(cell_color_vec3, status_level)
  DATA$group <- factor(DATA$status, levels = status_level)
  scores <- grep("auc_", colnames(DATA), value = TRUE)
  allPlotList <- list()
  
  if(plot_mode == "cell_facet"){
    comp_list <- list(c("Healthy controls","Polygenic JIA"),
                      c("Healthy controls","Monogenic JIA"),
                      c("Polygenic JIA","Monogenic JIA"))
    cells <- unique(as.character(DATA$celltype))
    names(cells) <- cells
    
    for (cell in cells) {
      cellD <- dplyr::filter(DATA, celltype == cell)
      onecell <- list()
      for (s in scores) {
        data1 <- cellD[, c("group", s)]
        colnames(data1) <- c("group", "score")
        data1$group <- droplevels(data1$group)
        if(all(is.na(data1$score))) next
        if(length(unique(na.omit(data1$group)))<2) next
        
        y_min  <- min(data1$score, na.rm = TRUE)
        y_max  <- max(data1$score, na.rm = TRUE)
        y_pad  <- (y_max - y_min) * 0.15
        ylim_lower <- y_min - y_pad
        ylim_upper <- y_max + y_pad
        sig_pos <- ylim_upper + (y_max - y_min) * 0.06
        
        p <- ggplot(data1, aes(x = group, y = score, fill = group)) +
          geom_violin(trim = FALSE, scale = "width", width = 0.32, size = 0.35, color = "black", position = position_dodge(width = 0.7)) +
          geom_boxplot(width = 0.35, color = "black", outlier.shape = NA, fill = NA, position = position_dodge(width = 0.7)) +
          stat_compare_means(
            method = "wilcox.test", label = "p.signif",
            comparisons = comp_list, tip.length = 0, step.increase = 0.04,
            hide.ns = FALSE, size = 2.5, y.position = sig_pos
          ) +
          labs(title = cell, x = "", y = s) +
          scale_fill_manual(values = color_map) +
          scale_x_discrete(expand = c(0.12,0)) +
          coord_cartesian(ylim = c(ylim_lower, ylim_upper), clip = "off") +
          theme_classic() +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            axis.text = element_text(colour = "black"),
            axis.text.x = element_blank(),
            axis.ticks.x=element_blank()
          )
        onecell[[s]] <- p
      }
      if(length(onecell)>0){
        allPlotList[[cell]] <- wrap_plots(onecell, ncol = length(onecell))
      }
    }
  } else if(plot_mode == "score_single"){
    names(scores) <- scores
    for(s in scores){
      data1 <- DATA[,c("celltype","group",s)]
      colnames(data1) <- c("celltype","group","score")
      
      p <- ggplot(data1, aes(x = celltype, y = score, fill = group)) +
        geom_violin(trim = FALSE, scale = "width", width = 0.32, size = 0.35, color = "black", position = position_dodge(width = 0.7)) +
        geom_boxplot(width = 0.35, color = "black", outlier.shape = NA, fill = NA, position = position_dodge(width = 0.7)) +
        labs(title = s, x = "", y = s) +
        scale_fill_manual(values = color_map) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      allPlotList[[s]] <- p
    }
  }
  return(allPlotList)
}

plt_cell <- Violin_Score_Adaptive_Sig(pbmc1, plot_mode = "cell_facet")
plt_cell$`CD14 Monocyte`
plt_cell$`CD14 Monocyte`/plt_cell$`CD16 Monocyte`/plt_cell$`Cytotoxic CD8 T`/plt_cell$pDC/plt_cell$`Memory B`
plt_cell$`Memory B`
plot_pdf_3cell_per_page <- function(plot_list, out_pdf){
  library(patchwork)
  cell_names <- names(plot_list)
  n_cell <- length(cell_names)
  pages <- ceiling(n_cell / 3)
  
  pdf(out_pdf, width = 14, height = 12)
  for(p in seq_len(pages)){
    idx_start <- (p-1)*3 + 1
    idx_end <- min(p*3, n_cell)
    pick <- cell_names[idx_start:idx_end]
    sub_plots <- plot_list[pick]
    # one column, each cell one row
    page_p <- wrap_plots(sub_plots, ncol = 1)
    print(page_p)
  }
  dev.off()
}
plt_cell <- Violin_Score_Adaptive_Sig(pbmc1, plot_mode = "cell_facet")
plot_pdf_3cell_per_page(plt_cell, out_pdf = "cell_violin_output.pdf")

###(3) specific mutation as one violin plot
Violin_All_row <- function(pbmc1, score, cell_order = c("CD14 Monocyte","CD16 Monocyte","pDC","Memory B"),
                           violin_width = 0.5, boxplot_width = 0.35){
  library(ggplot2); library(dplyr); library(patchwork)
  DATA <- pbmc1@meta.data
  if(!score %in% colnames(DATA)) stop("score not found in metadata")
  all_levels <- c("HC","Polygenic","157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  all_colors <- c("HC"="#6699CC","Polygenic"="#FEC260",
                  "157LACC1"="#B5BBE3",
                  "182LPIN2"="#E2A9C9","181NOD2"="#CBDAA9","183NOD2"="#B1DA99",
                  "200PSTPIP1"="#7B1FA2","190PSTPIP1"="#FFEFC1")
  DATA$group <- NA
  DATA$group[DATA$datasets %in% paste0("C", 1:6)] <- "HC"
  DATA$group[DATA$datasets %in% paste0("poly", 1:6)] <- "Polygenic"
  mono_ids <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  DATA$group[DATA$datasets %in% mono_ids] <- as.character(DATA$datasets[DATA$datasets %in% mono_ids])
  DATA <- DATA[!is.na(DATA$group) & DATA$celltype %in% cell_order, ]
  DATA$group <- factor(DATA$group, levels = all_levels)
  DATA$celltype <- factor(DATA$celltype, levels = cell_order)
  comp_pairs <- combn(all_levels, 2, simplify = FALSE)
  stats_list <- list(); kw_list <- list(); sig_cells <- c()
  for(ct in cell_order){
    ct_data <- DATA[DATA$celltype == ct, ]
    kw <- kruskal.test(ct_data[[score]] ~ ct_data$group)
    kw_list[[ct]] <- data.frame(celltype = ct, score = score, kw_p = kw$p.value, stringsAsFactors = FALSE)
    for(pair in comp_pairs){
      g1 <- ct_data[[score]][ct_data$group == pair[1]]
      g2 <- ct_data[[score]][ct_data$group == pair[2]]
      if(length(g1) < 5 | length(g2) < 5 | all(is.na(g1)) | all(is.na(g2))) next
      wt <- wilcox.test(g1, g2)
      stats_list[[length(stats_list)+1]] <- data.frame(
        celltype = ct, group1 = pair[1], group2 = pair[2],
        n1 = length(g1), n2 = length(g2),
        median1 = median(g1, na.rm = TRUE), median2 = median(g2, na.rm = TRUE),
        p_value = wt$p.value, stringsAsFactors = FALSE)
    }
  }
  kw_df <- do.call(rbind, kw_list)
  kw_df$kw_padj <- p.adjust(kw_df$kw_p, method = "BH")
  fmt_p <- function(p){
    if(p < 2.2e-16) return("< 2.2e-16")
    formatC(p, format = "e", digits = 2)
  }
  plot_list <- list()
  for(i in seq_along(cell_order)){
    ct <- cell_order[i]
    ct_data <- DATA[DATA$celltype == ct, ]
    kw_padj <- kw_df$kw_padj[kw_df$celltype == ct]
    p_str <- fmt_p(kw_padj)
    if(kw_padj < 0.05) sig_cells <- c(sig_cells, ct)
    hc_mean <- mean(ct_data[[score]][ct_data$group == "HC"], na.rm = TRUE)
    p <- ggplot(ct_data, aes(x = group, y = .data[[score]], fill = group)) +
      geom_violin(trim = FALSE, scale = "width", width = violin_width, color = NA) +
      geom_boxplot(width = boxplot_width, color = "black", outlier.shape = NA, fill = NA) +
      geom_hline(yintercept = hc_mean, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      scale_fill_manual(values = all_colors) +
      labs(title = bquote(atop(.(ct), italic(P) == .(p_str))),
           x = "", y = ifelse(i == 1, score, "")) +
      theme_classic() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9),
            plot.margin = margin(2, 2, 2, 2), legend.position = "none",
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.y = element_line(color = "black"),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.y = element_text(colour = "black", size = 9))
    plot_list[[ct]] <- p
  }
  if(length(stats_list) == 0){
    stats_df <- data.frame(celltype=character(), group1=character(), group2=character(),
                           n1=numeric(), n2=numeric(), median1=numeric(), median2=numeric(),
                           p_value=numeric(), p_adj=numeric(), stringsAsFactors = FALSE)
  } else {
    stats_df <- do.call(rbind, stats_list) %>% group_by(celltype) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>% arrange(celltype, p_adj) %>% as.data.frame()
  }
  combined <- wrap_plots(plot_list, ncol = length(cell_order)) + plot_layout(guides = "collect") &
    theme(plot.margin = margin(2, 2, 2, 2))
  cat("Score:", score, "| KW p.adj sig:", paste(sig_cells, collapse = ", "),
      "| Wilcoxon sig:", sum(stats_df$p_adj < 0.05), "\n")
  print(kw_df)
  return(list(plot = combined, kw = kw_df, stats = stats_df,
              sig = stats_df[stats_df$p_adj < 0.05, ], significant = sig_cells))
}

all_cells <- c("CD14 Monocyte","CD16 Monocyte","pDC","NK","Cytotoxic CD8 T","Naive B","Memory B")
res_nfkb <- Violin_All_row(pbmc1, score = "auc_NFKB", cell_order = all_cells, violin_width = 0.7, boxplot_width = 0.715)
res_mapk <- Violin_All_row(pbmc1, score = "auc_MAPK", cell_order = all_cells, violin_width = 0.7, boxplot_width = 0.715)
res_nfkb$plot / res_mapk$plot

all_cells2 <- c("CD14 Monocyte","CD16 Monocyte","pDC")
res_ifn2 <- Violin_All_row(pbmc1, score = "auc_IFNII", cell_order = all_cells2, violin_width = 0.7, boxplot_width = 0.715)
res_ifn1 <- Violin_All_row(pbmc1, score = "auc_IFNI", cell_order = all_cells2, violin_width = 0.7, boxplot_width = 0.715)
res_ifn1$plot | res_ifn2$plot
sig_nfkb <- res_nfkb$sig; sig_nfkb$score <- "auc_NFKB"
sig_mapk <- res_mapk$sig; sig_mapk$score <- "auc_MAPK"
sig_ifn2 <- res_ifn2$sig; sig_ifn2$score <- "auc_IFNII"
sig_ifn1 <- res_ifn1$sig; sig_ifn1$score <- "auc_IFNI"
all_sig <- rbind(sig_nfkb, sig_mapk, sig_ifn2,sig_ifn1)
all_sig <- all_sig[, c("score","celltype","group1","group2","n1","n2","median1","median2","p_value","p_adj")]
openxlsx::write.xlsx(all_sig, "E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\violin_wilcox_sig.xlsx", rowNames = FALSE)

###(4) per-sample mean AUCell score, boxplot by 3 groups, all group points, cell separators
Score_Group_Box <- function(pbmc1, score,
                            cell_levels = c("Naive CD4 T","Memory CD4 T","Naive CD8 T","Cytotoxic CD8 T","NK",
                                            "Treg","γδT","Naive B","Memory B","Plasma","CD14 Monocyte",
                                            "CD16 Monocyte","pDC"),
                            hc_col = "gray90", poly_col ="#F1BAAB", mono_col = "#D9C2D9",
                            hc_pt_col = "gray50", poly_pt_col = "#4E8FB8",
                            box_col = "black", box_lwd = 0.3, box_fatten = 1,
                            mono_colors = c("157LACC1"="#8B9CD8","182LPIN2"="#C47CA6","181NOD2"="#A5BE76",
                                            "183NOD2"="#83BE62","200PSTPIP1"="#5C1778","190PSTPIP1"="#EFD97E"),
                            mono_shapes = c("157LACC1"=16,"182LPIN2"=17,"181NOD2"=15,
                                            "183NOD2"=18,"200PSTPIP1"=8,"190PSTPIP1"=7),
                            box_width = 0.24, grp_off = 0.25, cell_gap = 0.3,
                            pt_size = 1.4, hc_poly_pt_size = 0.7, pt_alpha = 0.8,
                            jitter_w = NULL, vline_col = "gray70", seed = 1){
  library(ggplot2); library(dplyr); library(tidyr)
  if(!is.null(seed)) set.seed(seed)
  DATA <- pbmc1@meta.data
  if(!score %in% colnames(DATA)) stop("score not found in metadata")
  hc_ids <- paste0("C", 1:6); poly_ids <- paste0("poly", 1:6); mono_ids <- names(mono_colors)
  smry <- DATA %>% filter(!is.na(.data[[score]])) %>%
    group_by(datasets, celltype) %>%
    summarise(value = mean(.data[[score]], na.rm = TRUE), .groups = "drop")
  smry$group <- NA
  smry$group[smry$datasets %in% hc_ids] <- "HC"
  smry$group[smry$datasets %in% poly_ids] <- "Poly"
  smry$group[smry$datasets %in% mono_ids] <- "Mono"
  smry <- smry[!is.na(smry$group), ]
  smry$gene_label <- ifelse(smry$group == "Mono", smry$datasets, NA)
  smry$celltype <- factor(smry$celltype, levels = cell_levels)
  smry$group <- factor(smry$group, levels = c("HC","Poly","Mono"))
  if(is.null(jitter_w)) jitter_w <- box_width * 0.5
  ct_center <- (seq_along(cell_levels) - 1) * (1 + cell_gap) + 1
  smry$xcenter <- ct_center[as.integer(smry$celltype)]
  smry$xpos <- smry$xcenter + ifelse(smry$group == "HC", -grp_off,
                                     ifelse(smry$group == "Poly", 0, grp_off))
  get_sig <- function(d, g1, g2){
    v1 <- d$value[d$group == g1]; v2 <- d$value[d$group == g2]
    if(length(v1) < 2 || length(v2) < 2) return("ns")
    wt <- wilcox.test(v1, v2, exact = FALSE); p <- wt$p.value
    if(is.na(p)) return("ns")
    if(p < 0.001) return("***"); if(p < 0.01) return("**"); if(p < 0.05) return("*"); return("ns")
  }
  all_ct <- levels(smry$celltype)
  sig_rows <- do.call(rbind, lapply(all_ct, function(ct){
    d <- smry[smry$celltype == ct, ]
    data.frame(celltype = ct,
               HC_vs_Mono = get_sig(d, "HC","Mono"),
               HC_vs_Poly = get_sig(d, "HC","Poly"),
               Poly_vs_Mono = get_sig(d, "Poly","Mono"), stringsAsFactors = FALSE)
  }))
  sig_rows <- sig_rows[sig_rows$HC_vs_Mono != "ns" | sig_rows$HC_vs_Poly != "ns" | sig_rows$Poly_vs_Mono != "ns", ]
  p <- ggplot(smry, aes(x = xpos, y = value)) +
    geom_vline(xintercept = ct_center[-length(ct_center)] + (1 + cell_gap)/2,
               linetype = "dashed", color = vline_col, linewidth = 0.3) +
    geom_boxplot(aes(group = interaction(celltype, group), fill = group),
                 width = box_width, outlier.shape = NA, color = box_col,
                 linewidth = box_lwd, fatten = box_fatten) +
    geom_point(data = subset(smry, group == "HC"),
               color = hc_pt_col, shape = 16,
               position = position_jitter(width = jitter_w, height = 0),
               size = hc_poly_pt_size, alpha = pt_alpha) +
    geom_point(data = subset(smry, group == "Poly"),
               color = poly_pt_col, shape = 16,
               position = position_jitter(width = jitter_w, height = 0),
               size = hc_poly_pt_size, alpha = pt_alpha) +
    geom_point(data = subset(smry, group == "Mono"),
               aes(color = gene_label, shape = gene_label),
               position = position_jitter(width = jitter_w, height = 0),
               size = pt_size, alpha = pt_alpha) +
    scale_shape_manual(values = mono_shapes) +
    scale_color_manual(values = mono_colors) +
    scale_fill_manual(values = c("HC" = hc_col, "Poly" = poly_col, "Mono" = mono_col)) +
    scale_x_continuous(breaks = ct_center, labels = cell_levels, expand = c(0.02, 0.02)) +
    labs(x = "", y = score) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.title = element_blank())
  return(list(plot = p, sig = sig_rows))
}

res_NFKB <- Score_Group_Box(pbmc1, 
                            poly_col = "#F1BAAB",
                            poly_pt_col =  "#EC706E",
                            score = "auc_NFKB", cell_gap = 0.42,
                            box_width = 0.4, grp_off = 0.41, box_lwd = 0.3,
                            box_fatten = 1.2, pt_size = 1.5,
                            hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02, seed = 1)
print(res_NFKB$plot); print(res_NFKB$sig, row.names = FALSE)

res_Inflam <- Score_Group_Box(pbmc1, 
                              poly_col = "#F1BAAB",
                              poly_pt_col =  "#EC706E",
                              score = "auc_Inflammatory", cell_gap = 0.42,
                             box_width = 0.4, grp_off = 0.41, box_lwd = 0.3,
                             box_fatten = 1.2, pt_size = 1.5,
                             hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02, seed = 1)
print(res_Inflam$plot); print(res_Inflam$sig, row.names = FALSE)

##(a) main figure
res_NFKB$plot+theme(legend.position = "none")|res_Inflam$plot
(res_NFKB$plot+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()))/res_Inflam$plot

# res_IFNI <- Score_Group_Box(pbmc1, score = "auc_IFNI", cell_gap = 0.33,
#                             box_width = 0.27, grp_off = 0.28, box_lwd = 0.3,
#                             box_fatten = 1.2, pt_size = 1.5,
#                             hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02, seed = 1)
# print(res_IFNI$plot); print(res_IFNI$sig, row.names = FALSE)

res_IFNII <- Score_Group_Box(pbmc1, 
                             poly_col = "#F1BAAB",
                             poly_pt_col =  "#EC706E",
                             score = "auc_IFNII", cell_gap = 0.42,
                             box_width = 0.4, grp_off = 0.41, box_lwd = 0.3,
                             box_fatten = 1.2, pt_size = 1.5,
                             hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02, seed = 1)
print(res_IFNII$plot); print(res_IFNII$sig, row.names = FALSE)


res_MAPK <- Score_Group_Box(pbmc1, 
                            poly_col = "#F1BAAB",
                            poly_pt_col =  "#EC706E",
                            score = "auc_MAPK", cell_gap = 0.42,
                            box_width = 0.4, grp_off = 0.41, box_lwd = 0.3,
                            box_fatten = 1.2, pt_size = 1.5,
                            hc_poly_pt_size = 0.7, pt_alpha = 0.95, jitter_w = 0.02, seed = 1)
print(res_MAPK$plot); print(res_MAPK$sig, row.names = FALSE)

##(b) suplementary figure
(res_IFNII$plot+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()))/res_MAPK$plot

#####2 ----Immune signature score by PROGENY-----
library(progeny)
library(AUCell)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(cowplot)
library(RColorBrewer)
library(scales)

#ProgEscore <- progeny(pbmc1, scale = TRUE, organism = "Human")
#ProgEscore <- cbind(ProgEscore, pbmc1@meta.data[match(rownames(ProgEscore),rownames(pbmc1@meta.data)),c('celltype','sample')])
#saveRDS(ProgEscore,file="E:/Cohort PPT/JIA/code/CellPreprocess/Progeny.RDS")
#pathway_cols <- c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK", "NFkB", "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT")
ProgEscore<-readRDS("E:/Cohort PPT/JIA/code/CellPreprocess/Progeny.RDS")


####(1) pheatmap  no color bar
sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
col2 <- colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
bk <- unique(c(seq(-1.5, 1.5, length = 100)))
plot_pathway_heatmap <- function(data, cells, pathways, col_palette, bk_vals, cell_width = 10, cell_height = 12.5, sample_order = sample_level){
  allheat <- list()
  for(i in seq_along(cells)){
    subdata <- data %>% filter(celltype == cells[i])
    avg_data <- subdata %>% pivot_longer(cols = all_of(pathways), names_to = "Pathway", values_to = "Score") %>%
      group_by(sample, Pathway) %>% summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = sample, values_from = mean_score)
    mat <- avg_data %>% column_to_rownames("Pathway") %>% as.matrix()
    mat <- mat[, intersect(sample_order, colnames(mat)), drop = FALSE]
    allheat[[i]] <- as.ggplot(pheatmap(mat, scale = "row", cluster_cols = F, cluster_rows = T, show_rownames = T,
                                       border_color = NA, main = cells[i], fontsize = 6, legend = T, treeheight_col = 0, treeheight_row = 0,
                                       legend_breaks = seq(-2,2,0.5), color = col_palette, annotation_legend = T, breaks = bk_vals,
                                       cellwidth = cell_width, cellheight = cell_height, fontsize_col = 8, fontsize_row = 10,
                                       clustering_distance_rows = "euclidean", clustering_method = "median"))
  }
  names(allheat) <- cells
  return(allheat)
}
col2<-colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
bk <- unique(c(seq(-1.5,1.5, length=100)))
pathways_sel <- c("JAK-STAT","NFkB","MAPK","TNFa","Hypoxia")
cells_sel <- c("CD14 Monocyte","CD16 Monocyte","pDC")

Prog_heats <- plot_pathway_heatmap(data = ProgEscore, cells = cells_sel, 
                                   pathways = pathways_sel,cell_width = 8, cell_height = 10.5,
                                   col_palette = col2, bk_vals = bk, sample_order = sample_level)
plot_grid(plotlist = Prog_heats, ncol = 1, align = "hv")

####(2) pheatmap color bar
library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales)

sample_level <- c(paste0("C", 1:6), paste0("poly", 1:6),
                  "157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
mono_ids <- c("157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")

plot_pathway_heatmap_gg <- function(data, cells, pathways, sample_order = sample_level,
                                    pal = colorRampPalette(c("#542788","#F7F7F7","#D6604D"))(101),
                                    z_lim = 2.5, gap_width = 0.10,
                                    group_colors = c("HC"="gray90", poly_col ="#F1BAAB", mono_col = "#D9C2D9"),
                                    cluster_rows = TRUE, clustering_method = "ward.D2",
                                    scale_mode = "row", border_col = "white", border_lwd = 0.4,
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
    p <- ggplot(df, aes(x = xpos, y = ypos, fill = value)) +
      geom_tile(width = 1, height = 1, color = border_col, linewidth = border_lwd) +
      scale_fill_gradientn(colors = pal, limits = c(-z_lim, z_lim), oob = squish, name = legend_title) +
      scale_x_continuous(breaks = as.numeric(xp), labels = cols, expand = c(0, 0)) +
      scale_y_continuous(breaks = n_row:1, labels = row_ord, expand = c(0, 0)) +
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
            plot.margin = margin(4, 8, 4, 2))
    for(k in seq_along(cols)){
      p <- p + annotate("rect", xmin = xp[k] - 0.5, xmax = xp[k] + 0.5,
                        ymin = n_row + 0.5, ymax = n_row + 1.5,
                        fill = group_colors[grp[k]], color = NA)
    }
    out[[ct]] <- p
  }
  names(out) <- cells
  out
}

ProgEscore <- readRDS("E:/Cohort PPT/JIA/code/CellPreprocess/Progeny.RDS")
pathways_sel <- c("JAK-STAT","NFkB","MAPK","TNFa","Hypoxia")
cells_sel <- c("CD14 Monocyte","CD16 Monocyte","pDC")

Prog_heats <- plot_pathway_heatmap_gg(data = ProgEscore, cells = cells_sel, pathways = pathways_sel,
                                      z_lim = 2.5, gap_width = 0.10)
wrap_plots(Prog_heats, ncol = 1)


#####3 ----Inflmmation gene visualization by violin plots-----
###(1) one panel for single celltype, single gene
VioplotRun <- function(sce1, genes, cells = c("CD14 Monocyte","CD16 Monocyte","pDC"),
                       mode = c("seurattype","ggplottype"),
                       pdf_file = NULL, stats_file = NULL,
                       n_col = 4, n_row = 4, width = 14, height = 14,
                       assay = "RNA", sample_col = "sample",
                       group_levels = c("HC","Poly","157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1"),
                       group_colors = c("HC" = "#6699CC","Poly" = "#FEC260","157LACC1" = "#B5BBE3","182LPIN2" = "#E2A9C9","181NOD2" = "#CBDAA9","183NOD2" = "#B1DA99","200PSTPIP1" = "#7B1FA2","190PSTPIP1" = "#FFEFC1")){
  library(Seurat); library(ggplot2); library(patchwork)
  mode <- match.arg(mode)
  rn <- rownames(sce1@meta.data)
  md <- sce1@meta.data
  s <- as.character(md[[sample_col]])
  md$plot_sample <- factor(ifelse(s %in% paste0("C", 1:6), "HC", ifelse(grepl("^poly|^ploy", s), "Poly", s)), levels = group_levels)
  rownames(md) <- rn
  sce1@meta.data <- md
  genes <- genes[genes %in% rownames(sce1[[assay]])]
  cells <- intersect(cells, unique(as.character(sce1@meta.data$celltype)))
  if(length(genes) == 0 || length(cells) == 0) stop("no valid genes or cells")
  expr <- tryCatch(GetAssayData(sce1, assay = assay, layer = "data"), error = function(e) GetAssayData(sce1, assay = assay, slot = "data"))
  sub_list <- lapply(cells, function(ct) subset(sce1, celltype == ct))
  names(sub_list) <- cells
  Vioplot <- list(); stats_rows <- list()
  for(ct in cells){
    sub_sce <- sub_list[[ct]]
    if(ncol(sub_sce) < 3) next
    ps <- as.character(sub_sce@meta.data$plot_sample)
    for(g in genes){
      key <- paste0(g, "|", ct)
      ttl <- paste(ct, g, sep = " | ")
      v <- as.numeric(expr[g, colnames(sub_sce)])
      if(mode == "seurattype"){
        p <- VlnPlot(sub_sce, features = g, group.by = "plot_sample", pt.size = 0, cols = group_colors) +
          theme_classic() + xlab(NULL) + ylab(NULL) + ggtitle(ttl) +
          theme(legend.position = "none", plot.title = element_text(size = 9, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
                axis.text.y = element_text(size = 10, colour = "black"),
                axis.ticks.x = element_blank(), plot.margin = margin(2, 4, 2, 2))
      } else {
        V1 <- data.frame(plot_sample = factor(ps, levels = group_levels), value = v)
        p <- ggplot(V1, aes(x = plot_sample, y = value)) +
          geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = plot_sample)) +
          scale_fill_manual(values = group_colors) + theme_classic() +
          xlab(NULL) + ylab(NULL) + ggtitle(ttl) +
          theme(legend.position = "none", plot.title = element_text(size = 9, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
                axis.text.y = element_text(size = 10, colour = "black"),
                axis.ticks.x = element_blank(),
                panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
                plot.margin = margin(4, 4, 4, 4))
      }
      Vioplot[[key]] <- p
      for(gp in group_levels){
        vv <- v[ps == gp]
        if(length(vv) < 2) next
        stats_rows[[length(stats_rows) + 1]] <- data.frame(gene = g, celltype = ct, group = gp,
                                                           n = length(vv), mean_expr = mean(vv), median_expr = median(vv),
                                                           pct_expressing = mean(vv > 0) * 100, stringsAsFactors = FALSE)
      }
    }
  }
  stats <- do.call(rbind, stats_rows)
  if(!is.null(stats_file) && nrow(stats) > 0) write.csv(stats, stats_file, row.names = FALSE)
  if(!is.null(pdf_file) && length(Vioplot) > 0){
    per_page <- n_col * n_row
    pdf(pdf_file, width = width, height = height)
    grp_idx <- split(seq_along(Vioplot), ceiling(seq_along(Vioplot) / per_page))
    for(i in grp_idx) print(wrap_plots(Vioplot[i], ncol = n_col, nrow = n_row))
    dev.off()
    cat("pdf:", pdf_file, "| plots:", length(Vioplot), "| pages:", length(grp_idx), "| stats:", nrow(stats), "\n")
  }
  Vioplot
}
load("E:\\通路基因集合\\gene4pathway.Rdata")
IFN28 <- c("DDX60","EPSTI1","HERC5","HERC6","IFI27","IFI44","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFIT5","ISG15","LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SIGLEC1","SPATS2L","USP18","CXCL10","GBP1","SOCS1")
pyrop <- c("AIM2","APIP","CASP1","CASP4","CASP8","DHX9","ELANE","GSDMA","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB","NAIP","NLRC4","NLRP1","NLRP6","NLRP9","TREM2","ZBP1","NLRP3","PYCARD","NEK7","IL1B")
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
auto_genes <- reactome %>% filter(gs_name == "REACTOME_AUTOPHAGY") %>% pull(gene_symbol) %>% unique()
hm <- msigdbr(species = "Homo sapiens", category = "H")
gly_genes <- hm %>% filter(gs_name == "HALLMARK_GLYCOLYSIS") %>% pull(gene_symbol) %>% unique()
oxphos_genes <- hm %>% filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% pull(gene_symbol) %>% unique()
arachidonic_genes <- c("PLA2G4A","PLA2G6","PLA2G2A","PTGS1","PTGS2","ALOX5","ALOX5AP","ALOX12","ALOX15","ALOX15B","LTC4S","PTGES","PTGES2","TBXAS1","GPX4","CYP2J2","CYP4F2","EPHX2","HPGDS")
inflam_genes <- hm %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique()


sets <- list(IFN28 = IFN28, pyrop = pyrop, autophagy = auto_genes, glycolysis = gly_genes,
             OXPHOS = oxphos_genes, arachidonic = arachidonic_genes, inflammatory = inflam_genes)

cells_use <- c("CD14 Monocyte","CD16 Monocyte","pDC")
out_dir <- "E:/Cohort PPT/JIA/code/CellPreprocess"
mode_lab <- c("seurattype" = "seu", "ggplottype" = "gg")

vio_res <- list()
for(md in names(mode_lab)){
  for(nm in names(sets)){
    key <- paste0(nm, "_", mode_lab[md])
    vio_res[[key]] <- VioplotRun(pbmc1, genes = sets[[nm]], cells = cells_use, mode = md,
                                 pdf_file = file.path(out_dir, paste0("violin_", nm, "_3cell_", mode_lab[md], ".pdf")),
                                 stats_file = file.path(out_dir, paste0("violin_", nm, "_3cell_", mode_lab[md], "_stats.csv")))
  }
}
names(vio_res)


###(2) one panel, single gene and multiple cell types
Vln_Seurat_AllCell <- function(seu_obj, genes, celltypes = NULL, assay = "RNA", pdf_file = NULL,
                               per_page = 4, width = 16, height = 14,
                               sample_col = "sample",
                               group_levels = c("HC","Poly","157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1"),
                               group_colors = c("HC" = "#6699CC","Poly" = "#FEC260","157LACC1" = "#B5BBE3","182LPIN2" = "#E2A9C9","181NOD2" = "#CBDAA9","183NOD2" = "#B1DA99","200PSTPIP1" = "#7B1FA2","190PSTPIP1" = "#FFEFC1")){
  library(Seurat); library(ggplot2); library(patchwork)
  celltype_levels <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T","γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC","Platelet","Erythrocyte")
  md <- seu_obj@meta.data
  s <- as.character(md[[sample_col]])
  md$plot_sample <- factor(ifelse(s %in% paste0("C", 1:6), "HC", ifelse(grepl("^poly|^ploy", s), "Poly", s)), levels = group_levels)
  md$celltype <- factor(as.character(md$celltype), levels = celltype_levels)
  rownames(md) <- rownames(seu_obj@meta.data)
  seu_obj@meta.data <- md
  if(!is.null(celltypes)){
    ct_all <- unique(as.character(seu_obj@meta.data$celltype))
    ct_use <- intersect(celltypes, ct_all)
    if(length(ct_use) == 0) stop("No valid celltypes")
    if(length(ct_use) < length(ct_all)) seu_obj <- subset(seu_obj, celltype %in% ct_use)
    seu_obj@meta.data$celltype <- factor(as.character(seu_obj@meta.data$celltype), levels = ct_use)
  }
  genes <- genes[genes %in% rownames(seu_obj[[assay]])]
  if(length(genes) == 0) stop("No valid genes in assay")
  plot_list <- list()
  for(g in genes){
    p <- VlnPlot(seu_obj, features = g, group.by = "celltype", split.by = "plot_sample", assay = assay, pt.size = 0) +
      scale_fill_manual(values = group_colors) +
      labs(title = g, x = "", y = g) + theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = "black"))
    plot_list[[g]] <- p
  }
  if(!is.null(pdf_file) && length(plot_list) > 0){
    pdf(pdf_file, width = width, height = height)
    grp_idx <- split(seq_along(plot_list), ceiling(seq_along(plot_list) / per_page))
    for(i in grp_idx) print(wrap_plots(plot_list[i], ncol = 1))
    dev.off()
    cat("pdf:", pdf_file, "| genes:", length(plot_list), "| pages:", length(grp_idx), "\n")
  }
  plot_list
}

pl_3cell <- Vln_Seurat_AllCell(pbmc1, genes = c(
  "NLRP3","NEK7","IFI44","IFI44L","DDX58","OAS3"),
  celltypes = c("CD14 Monocyte","CD16 Monocyte","pDC"),
  per_page = 2, width = 11, height = 9)
library(patchwork)
fig <- wrap_plots(pl_3cell, ncol = 2, byrow = TRUE)
fig


