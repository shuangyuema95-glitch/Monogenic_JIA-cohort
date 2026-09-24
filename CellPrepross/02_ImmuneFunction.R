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
pyrop<-c("AIM2","APIP","CASP1","CASP4","CASP8","DHX9","ELANE","GSDMA","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB","NAIP","NLRC4","NLRP1","NLRP6","NLRP9","TREM2","ZBP1")#form Spectra
reactome <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:REACTOME")
IL1_genes <- reactome %>%filter(gs_name == "REACTOME_INTERLEUKIN_1_SIGNALING") %>%pull(gene_symbol) %>%unique()

pbmc1 <- score_Calu(pbmc1, NFKB$gene, "NFKB")
pbmc1 <- score_Calu(pbmc1, IFN28, "IFNI")
pbmc1 <- score_Calu(pbmc1, IFN2, "IFNII")
pbmc1 <- score_Calu(pbmc1, MAPK$gene, "MAPK")
pbmc1 <- score_Calu(pbmc1, pyrop, "pyrop")
pbmc1 <- score_Calu(pbmc1,IL1_genes, "IL1")
scores<-pbmc1@meta.data;saveRDS(scores,file="0827_pbmc1_aucScore.rds")
##scores loading
scores<-readRDS("0827_pbmc1_aucScore.rds")
pbmc1@meta.data<-scores

###(2)boxplot visualization
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
          geom_violin(trim = FALSE, scale = "width", width = 0.32, size = 0.35, color = "black") +
          geom_boxplot(width = 0.35, color = "black", outlier.shape = NA, fill = NA) +
          stat_compare_means(
            method = "wilcox.test", label = "p.signif",
            comparisons = comp_list, tip.length = 0, step.increase = 0.04,
            hide.ns = FALSE, size = 2.5, y.position = sig_pos
          ) +
          labs(title = cell, x = "", y = s) +
          scale_fill_manual(values = color_map) +
          coord_cartesian(ylim = c(ylim_lower, ylim_upper), clip = "off") +
          theme_classic() +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            #axis.ticks = element_line(color = "black"),
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
    ## one score,one plot; x‑axis = celltype, fill=group(status), NO significance
    names(scores) <- scores
    for(s in scores){
      data1 <- DATA[,c("celltype","group",s)]
      colnames(data1) <- c("celltype","group","score")
      
      p <- ggplot(data1, aes(x = celltype, y = score, fill = group)) +
        geom_violin(trim = FALSE, scale = "width", width = 0.32, size = 0.35, color = "black") +
        geom_boxplot(width = 0.35, color = "black", outlier.shape = NA, fill = NA) +
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
# plt_cell$NK
plt_cell$`Memory B`
# plt_cell$`Naive B`
#plt_cell$Treg
# plt_cell$pDC
#plt_cell$Plasma
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

###(3)specific mutation as one violin plot
###### Violin_All_row
###### Violin_All_row
###### Violin_All_row
Violin_All_row <- function(pbmc1, score, cell_order = c("CD14 Monocyte","CD16 Monocyte","pDC","Memory B"),
                           violin_width = 0.5, boxplot_width = 0.35){
  library(ggplot2); library(dplyr); library(patchwork)
  DATA <- pbmc1@meta.data
  if(!score %in% colnames(DATA)) stop("score not found in metadata")
  all_levels <- c("HC","Polygenic","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  all_colors <- c("HC"="#6699CC","Polygenic"="#FEC260",
                  "182LPIN2"="#E2A9C9","181NOD2"="#CBDAA9","183NOD2"="#B1DA99",
                  "200PSTPIP1"="#7B1FA2","190PSTPIP1"="#FFEFC1")
  DATA$group <- NA
  DATA$group[DATA$datasets %in% paste0("C",1:6)] <- "HC"
  DATA$group[DATA$datasets %in% paste0("ploy",1:6)] <- "Polygenic"
  mono_ids <- c("181NOD2","182LPIN2","183NOD2","190PSTPIP1","200PSTPIP1")
  DATA$group[DATA$datasets %in% mono_ids] <- DATA$datasets[DATA$datasets %in% mono_ids]
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
res_nfkb$plot / res_mapk$plot #main figure

all_cells<-c("CD14 Monocyte","CD16 Monocyte","pDC")
res_ifn2 <- Violin_All_row(pbmc1, score = "auc_IFNII", cell_order = all_cells, violin_width = 0.7, boxplot_width = 0.715)
res_ifn1 <- Violin_All_row(pbmc1, score = "auc_IFNI", cell_order = all_cells, violin_width = 0.7, boxplot_width = 0.715)
res_ifn1$plot | res_ifn2$plot #supplementary figure

sig_nfkb <- res_nfkb$sig; sig_nfkb$score <- "auc_NFKB"
sig_mapk <- res_mapk$sig; sig_mapk$score <- "auc_MAPK"
sig_ifn2 <- res_ifn2$sig; sig_ifn2$score <- "auc_IFNII"
sig_ifn1 <- res_ifn1$sig; sig_ifn1$score <- "auc_IFNI"
all_sig <- rbind(sig_nfkb, sig_mapk, sig_ifn2,sig_ifn1)
all_sig <- all_sig[, c("score","celltype","group1","group2","n1","n2","median1","median2","p_value","p_adj")]
openxlsx::write.xlsx(all_sig, "E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\violin_wilcox_sig.xlsx", rowNames = FALSE)

#####2 ----Immune signature score by PROGENY-----
library(progeny)
library(AUCell)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(cowplot)
library(RColorBrewer)
library(scales)

sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")

plot_pathway_heatmap <- function(data, cells, col_palette, bk_vals){
  pathway_cols <- c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK", "NFkB", "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT")
  allheat <- list()
  for(i in seq_along(cells)){
    subdata <- data %>% filter(celltype == cells[i])
    avg_data <- subdata %>% pivot_longer(cols = all_of(pathway_cols), names_to = "Pathway", values_to = "Score") %>%
      group_by(sample, Pathway) %>% summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = sample, values_from = mean_score)
    mat <- avg_data %>% column_to_rownames("Pathway") %>% as.matrix()
    mat <- mat[,intersect(sample_level,colnames(mat)),drop=FALSE]
    allheat[[i]] <- as.ggplot(pheatmap(mat, scale = "row", cluster_cols = F, cluster_rows = T, show_rownames = T,
                                       border_color = NA, main = cells[i], fontsize = 6, legend = T, treeheight_col = 0, treeheight_row = 0,
                                       legend_breaks = seq(-2,2,0.5), color = col_palette, annotation_legend = T, breaks = bk_vals,
                                       cellwidth = 10, cellheight = 12.5, fontsize_col = 8, fontsize_row = 10,
                                       clustering_distance_rows = "euclidean", clustering_method = "median"))
  }
  return(allheat)
}

#ProgEscore <- progeny(pbmc1, scale = TRUE, organism = "Human")
#ProgEscore <- cbind(ProgEscore, pbmc1@meta.data[match(rownames(ProgEscore),rownames(pbmc1@meta.data)),c('celltype','sample')])
#saveRDS(ProgEscore,file="E:/Cohort PPT/JIA/code/CellPreprocess/Progeny.RDS")
col2<-colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
#col2<-heat_col <- colorRampPalette(c("#6699CC", "#FFFFFF", "#AF478A"))(101)
#col2<-heat_col <- colorRampPalette(c("#007ABA", "#FFFFFF", "#BB0021"))(101)
#col2<- colorRampPalette(c("#847AB3", "#FFFFFF", "#EC706E"))(101)

bk <- unique(c(seq(-1.5,1.5, length=100)))
cells <- unique(as.character(pbmc1$celltype))
Prog_heats <- plot_pathway_heatmap(ProgEscore, cells, col_palette=col2, bk_vals=bk)
plot_grid(plotlist = Prog_heats, ncol = 3, align = "hv")


#####3 ----Inflmmation gene visualization by violin plots-----
###(1) total cells
Vln_Seurat_AllCell <- function(seu_obj, genes, assay = "RNA"){
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  
  sample_level <- c("C1","C2","C3","C4","C5","C6",
                    "poly1","poly2","poly3","poly4","poly5","poly6",
                    "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  
  celltype_levels <- c("Naive CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                       "γδT","Naive B","Memory B","Plasma","CD14 Monocyte",
                       "CD16 Monocyte","pDC","Platelet","Erythrocyte")
  
  seu_obj@meta.data <- seu_obj@meta.data %>%
    mutate(
      plot_sample = case_when(
        sample %in% c("C1","C2","C3","C4","C5","C6") ~ "HC",
        grepl("^poly", sample) ~ "Poly",
        TRUE ~ as.character(sample)
      ),
      celltype = factor(celltype, levels = celltype_levels)
    )
  
  final_group_levels <- c("HC","Poly","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  seu_obj@meta.data$plot_sample <- factor(seu_obj@meta.data$plot_sample, levels = final_group_levels)
  
  color_map <- c(
    "HC"        = "#6699CC",
    "Poly"      = "#FEC260",
    "182LPIN2"  = "#E2A9C9",
    "181NOD2"   = "#CBDAA9",
    "183NOD2"   = "#B1DA99",
    "200PSTPIP1"= "#7B1FA2",
    "190PSTPIP1"= "#FFEFC1"
  )
  
  genes <- genes[genes %in% rownames(GetAssay(seu_obj, assay = assay))]
  if(length(genes) == 0) stop("No valid genes found in given assay")
  
  plot_list <- list()
  names(genes) <- genes
  for(g in genes){
    p <- VlnPlot(
      seu_obj,
      features = g,
      group.by = "celltype",
      split.by = "plot_sample",
      assay = assay,
      pt.size = 0
    ) +
      scale_fill_manual(values = color_map) +
      labs(title = g, x = "", y = g) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black")
      )
    plot_list[[g]] <- p
  }
  return(plot_list)
}
plist <- Vln_Seurat_AllCell(pbmc1, genes = IFN28)
plist2 <- Vln_Seurat_AllCell(pbmc1, genes = pyrop)

###(2) single celltype, single gene
all_target_gene <- list(
  "CD14 Monocyte"    = c("EPSTI1", "IFI44"),
  "Cytotoxic CD8 T"  = c("EPSTI1", "NLRP1"),
  "NK"               = c("EPSTI1", "IFI6", "NLRP1"),
  "Treg"             = c("EPSTI1", "NLRP1"),
  "CD16 Monocyte"    = c("IFI44", "IFI44L", "IFIT2", "IFIT3", "NLRC4"),
  "pDC"              = c("IFI44", "IFI6"),
  "Plasma"           = c("GSDMB"))
target_gene_full <- list(
  list(gene = "EPSTI1",  cell = c("CD14 Monocyte","Cytotoxic CD8 T","NK","Treg","Memory B")),
  list(gene = "IFI44",   cell = c("CD14 Monocyte","CD16 Monocyte","pDC")),
  list(gene = "IFI6",    cell = c("NK","pDC")),
  list(gene = "IFI44L",  cell = c("CD16 Monocyte")),
  list(gene = "IFIT2",   cell = c("CD16 Monocyte")),
  list(gene = "IFIT3",   cell = c("CD16 Monocyte")),
  list(gene = "GSDMB",   cell = c("Plasma")),
  list(gene = "NLRC4",   cell = c("CD16 Monocyte")),
  list(gene = "NLRP1",   cell = c("NK","Cytotoxic CD8 T","Treg"))
)
VioplotRun <- function(sce1, mode = c("seurattype","ggplottype"), target_list){
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  
  mode <- match.arg(mode)
  
  sample_level <- c("C1","C2","C3","C4","C5","C6",
                    "poly1","poly2","poly3","poly4","poly5","poly6",
                    "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  
  final_groups <- c("HC","Poly","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
  
  color_vec <- c(
    "HC"         = "#6699CC",
    "Poly"       = "#FEC260",
    "182LPIN2"   = "#E2A9C9",
    "181NOD2"    = "#CBDAA9",
    "183NOD2"    = "#B1DA99",
    "200PSTPIP1" = "#7B1FA2",
    "190PSTPIP1" = "#FFEFC1"
  )
  
  sce1@meta.data <- sce1@meta.data %>%
    mutate(
      plot_sample = case_when(
        sample %in% c("C1","C2","C3","C4","C5","C6") ~ "HC",
        grepl("^poly", sample) ~ "Poly",
        TRUE ~ as.character(sample)
      ),
      plot_sample = factor(plot_sample, levels = final_groups)
    )
  
  Vioplot <- list()
  
  for(entry in target_list){
    g <- entry$gene
    cell_vec <- entry$cell
    
    for(ct in cell_vec){
      plot_key <- paste0(g,"|",ct)
      
      if(mode == "seurattype"){
        sub_sce <- subset(sce1, celltype == ct)
        if(ncol(sub_sce) < 3) next
        p <- VlnPlot(sub_sce, features = g, group.by = "plot_sample", pt.size = 0, cols = color_vec) +
          xlab(NULL) + ylab(NULL) + ggtitle(g) +
          theme(
            plot.margin = margin(0,0,0,0),
            legend.position = "none",
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 11),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 12, hjust = 0.5, lineheight = 0.2)
          )
        Vioplot[[plot_key]] <- p
        
      } else if(mode == "ggplottype"){
        data_exp <- GetAssayData(sce1, assay = "RNA", slot = "data")
        V1 <- sce1@meta.data %>% dplyr::filter(celltype == ct)
        if(nrow(V1) < 3) next
        idx_g <- match(g, rownames(data_exp))
        if(is.na(idx_g)) next
        V1$value <- data_exp[idx_g, match(rownames(V1), colnames(data_exp))]
        
        p <- ggplot(data = V1, mapping = aes(x = plot_sample, y = value)) +
          geom_violin(scale = "width", adjust = 1, trim = TRUE, mapping = aes(fill = plot_sample)) +
          scale_fill_manual(values = color_vec) +
          theme_classic() +
          xlab(NULL) + ylab(NULL) + ggtitle(paste(ct,g,sep = " ")) +
          theme(
            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
            legend.position = "none",
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.ticks.x = element_blank(),
            panel.border = element_rect(color = "black", size = 0.5, fill = NA),
            plot.title = element_text(size = 12, hjust = 0.5, lineheight = 0.2),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        Vioplot[[plot_key]] <- p
      }
    }
  }
  return(Vioplot)
}

target_gene_full <- list(
  list(gene = "EPSTI1",  cell = c("CD14 Monocyte","Cytotoxic CD8 T","Treg","Memory B","NK")),
  list(gene = "IFI44",   cell = c("CD14 Monocyte","CD16 Monocyte","pDC")),
  list(gene = "IFI6",    cell = c("pDC","NK")),
  list(gene = "IFI44L",  cell = c("CD16 Monocyte")),
  list(gene = "IFIT2",   cell = c("CD16 Monocyte")),
  list(gene = "IFIT3",   cell = c("CD16 Monocyte")),
  list(gene = "GSDMB",   cell = c("Plasma")),
  list(gene = "NLRC4",   cell = c("CD16 Monocyte")),
  list(gene = "NLRP1",   cell = c("Cytotoxic CD8 T","Treg","NK"))
)
target_gene_full2<-list(
  list(gene="TNF",cell=c("CD14 Monocyte","CD16 Monocyte")),
  list(gene="IL6",cell=c("CD14 Monocyte","CD16 Monocyte")),
  list(gene="IL18",cell=c("CD14 Monocyte","CD16 Monocyte")),
  list(gene="CXCL8",cell=c("CD14 Monocyte","CD16 Monocyte")),
  list(gene="IL1B",cell=c("CD14 Monocyte","CD16 Monocyte"))
  
)

vio_seu <- VioplotRun(sce1 = pbmc1, mode = "seurattype", target_list = target_gene_full2)
#vio_gg  <- VioplotRun(sce1 = pbmc1, mode = "ggplottype", target_list = target_gene_full)
vio_seu_out <- list()
for(k in names(vio_seu)){
  p <- vio_seu[[k]]
  parts <- strsplit(k,"\\|")[[1]]
  gnm <- parts[1]
  ctnm <- parts[2]
  p_new <- p +
    ggtitle(paste0(ctnm," | ",gnm)) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 8, face = "plain", hjust = 0.5))
  vio_seu_out[[k]] <- p_new
}
fig_all <- wrap_plots(vio_seu_out, ncol =6, nrow = 3)
fig_all

#####4 ----Differentially expressed genes across cell types-----
###(1)DEG analysis
library(tibble)
#cell_types <- c("Cytotoxic CD8 T","Memory B","CD14 Monocyte","CD16 Monocyte","pDC","NK","Plasma","Treg")
cell_types<-unique(pbmc1$celltype)
comp_list <- list(
  Poly_vs_HC = c("Polygenic JIA","Healthy controls"),
  Mono_vs_Poly = c("Monogenic JIA","Polygenic JIA"),
  Mono_vs_HC = c("Monogenic JIA","Healthy controls")
)

de_res_list <- list()
for(ct in cell_types){
  sub_pbmc1 <- subset(pbmc1, celltype == ct)
  if(ncol(sub_pbmc1) < 10){ de_res_list[[ct]] <- NULL; next }
  Idents(sub_pbmc1) <- "status"
  
  tmp_l <- list()
  for(nm in names(comp_list)){
    id1 <- comp_list[[nm]][1]
    id2 <- comp_list[[nm]][2]
    n1 <- sum(Idents(sub_pbmc1) == id1)
    n2 <- sum(Idents(sub_pbmc1) == id2)
    
    if(n1 < 3 || n2 < 3){
      message(sprintf("skip %s | %s: n(%s)=%d, n(%s)=%d", ct, nm, id1, n1, id2, n2))
      next
    }
    
    mk <- FindMarkers(
      sub_pbmc1,
      ident.1 = id1, ident.2 = id2,
      only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25
    )
    mk <- mk %>%
      rownames_to_column("gene") %>%
      mutate(comparison = nm) %>%
      filter(abs(avg_log2FC) > 1, p_val_adj < 0.01)
    tmp_l[[nm]] <- mk
  }
  
  if(length(tmp_l) > 0){
    de_res_list[[ct]] <- bind_rows(tmp_l)
  } else {
    de_res_list[[ct]] <- NULL
  }
}


library(dplyr)
gene_sets <- list(
  IFN28 = IFN28,
  IFN2 = IFN2$gene,
  MAPK = MAPK$gene,
  NFKB = NFKB$gene,
  pyrop = pyrop,
  IL1 = IL1_genes
)

de_res_filter <- de_res_list[!sapply(de_res_list,is.null)]
all_de_df <- NULL
for(ct_name in names(de_res_filter)){
  df_tmp <- de_res_filter[[ct_name]]
  df_tmp$celltype <- ct_name
  all_de_df <- rbind(all_de_df,df_tmp)
}
all_de_df <- as_tibble(all_de_df)
all_de_df <- all_de_df %>% mutate(direction=ifelse(avg_log2FC>0,"up","down"))

stat_df <- all_de_df %>%
  group_by(celltype,comparison,direction) %>%
  summarise(gene_count=n(),.groups="drop")

split_list <- split(all_de_df,paste(all_de_df$celltype,all_de_df$comparison,all_de_df$direction,sep="_"))
split_list <- split_list[sapply(split_list,nrow)>0]

intersect_res <- lapply(split_list,function(df){
  gvec <- df$gene
  out <- lapply(gene_sets,function(s) length(intersect(gvec,s)))
  tibble(!!!out) %>%
    mutate(celltype=unique(df$celltype),comparison=unique(df$comparison),direction=unique(df$direction))
}) %>% bind_rows()

summary_final <- left_join(stat_df,intersect_res,by=c("celltype","comparison","direction"))
View(summary_final)
write.csv(summary_final,"./CellPreprocess/de_gene_set_summary.csv",row.names=F)

summary_final2<-summary_final[summary_final$celltype%in%c("Naive B","Memory B",
"Naive CD4 T","Naive CD8 T","Cytotoxic CD8 T","γδT","Plasma","CD14 Monocyte","CD16 Monocyte"),]

library(tidyverse)
library(scatterpie)
plot_deg_pie_grid_summary <- function(
    dat,
    celltype_order,
    comparison_order,
    fixed_r = 0.8,
    col_up = "#D1929B",
    col_down = "#C2D4F2"
){
  
  deg_summary <- dat %>%
    pivot_wider(names_from = direction, values_from = gene_count, values_fill = 0) %>%
    group_by(celltype, comparison) %>%
    summarise(
      up_count = sum(up),
      down_count = sum(down),
      total_deg_ct = up_count + down_count,
      .groups = "drop"
    ) %>%
    mutate(
      celltype  = factor(celltype, levels = celltype_order),
      comparison = factor(comparison, levels = comparison_order)
    )
  
  # plot1: x=celltype, y=comparison
  coord_base1 <- deg_summary %>%
    mutate(
      x = as.integer(celltype)*3.2,
      y = as.integer(comparison)*3.2,
      r_size = fixed_r
    )
  coord_label1 <- coord_base1 %>%
    mutate(
      frac_up = up_count / total_deg_ct,
      frac_down = down_count / total_deg_ct,
      angle_up = frac_up * pi,
      lab_x_up = x + r_size * 0.5 * cos(angle_up),
      lab_y_up = y + r_size * 0.5 * sin(angle_up),
      angle_down = frac_up * pi + frac_down * pi,
      lab_x_down = x + r_size * 0.5 * cos(angle_down),
      lab_y_down = y + r_size * 0.5 * sin(angle_down)
    )
  full_x_breaks1 <- seq_along(celltype_order)*3.2
  full_y_breaks1 <- seq_along(comparison_order)*3.2
  
  p1 <- ggplot() +
    geom_scatterpie(data=coord_base1,aes(x=x,y=y,r=r_size),cols=c("up_count","down_count"),size=0,color=NA)+
    scale_fill_manual(values=c(col_up, col_down),labels=c("UP","DOWN"),name="DEG direction")+
    geom_text(data=filter(coord_label1, up_count>0), aes(x=lab_x_up,y=lab_y_up,label=up_count), size=2.4, color="black") +
    geom_text(data=filter(coord_label1, down_count>0), aes(x=lab_x_down,y=lab_y_down,label=down_count), size=2.4, color="black") +
    scale_x_continuous(breaks = full_x_breaks1,labels = celltype_order,expand = c(0,0))+
    scale_y_continuous(breaks = full_y_breaks1,labels = comparison_order,expand = c(0,0),trans = "reverse")+
    labs(x=NULL,y=NULL)+
    theme_bw()+
    theme(
      panel.grid.major.x = element_line(colour="grey75",linewidth=0.25,linetype="dashed"),
      panel.grid.major.y = element_line(colour="grey75",linewidth=0.25,linetype="22"),
      axis.text.x = element_text(angle=45,hjust=1,size=7.5,color="black"),
      axis.text.y = element_text(size=7,color="black"),
      axis.ticks = element_line(color="black"),
      axis.line = element_line(color="black"),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(2,2,2,2,"mm")
    )
  
  # plot2: swap x/y: x=comparison, y=celltype
  coord_base2 <- deg_summary %>%
    mutate(
      x = as.integer(comparison)*3.2,
      y = as.integer(celltype)*3.2,
      r_size = fixed_r
    )
  coord_label2 <- coord_base2 %>%
    mutate(
      frac_up = up_count / total_deg_ct,
      frac_down = down_count / total_deg_ct,
      angle_up = frac_up * pi,
      lab_x_up = x + r_size * 0.5 * cos(angle_up),
      lab_y_up = y + r_size * 0.5 * sin(angle_up),
      angle_down = frac_up * pi + frac_down * pi,
      lab_x_down = x + r_size * 0.5 * cos(angle_down),
      lab_y_down = y + r_size * 0.5 * sin(angle_down)
    )
  full_x_breaks2 <- seq_along(comparison_order)*3.2
  full_y_breaks2 <- seq_along(celltype_order)*3.2
  
  p2 <- ggplot() +
    geom_scatterpie(data=coord_base2,aes(x=x,y=y,r=r_size),cols=c("up_count","down_count"),size=0,color=NA)+
    scale_fill_manual(values=c(col_up, col_down),labels=c("UP","DOWN"),name="DEG direction")+
    geom_text(data=filter(coord_label2, up_count>0), aes(x=lab_x_up,y=lab_y_up,label=up_count), size=2.4, color="black") +
    geom_text(data=filter(coord_label2, down_count>0), aes(x=lab_x_down,y=lab_y_down,label=down_count), size=2.4, color="black") +
    scale_x_continuous(breaks = full_x_breaks2,labels = comparison_order,expand = c(0,0))+
    scale_y_continuous(breaks = full_y_breaks2,labels = celltype_order,expand = c(0,0),trans = "reverse")+
    labs(x=NULL,y=NULL)+
    theme_bw()+
    theme(
      panel.grid.major.x = element_line(colour="grey75",linewidth=0.25,linetype="dashed"),
      panel.grid.major.y = element_line(colour="grey75",linewidth=0.25,linetype="22"),
      axis.text.x = element_text(angle=45,hjust=1,size=7.5,color="black"),
      axis.text.y = element_text(size=7,color="black"),
      axis.ticks = element_line(color="black"),
      axis.line = element_line(color="black"),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(2,2,2,2,"mm")
    )
  
  return(list(p_celltype_x=p1, p_comparison_x=p2))
}

# usage
my_celltype_order <- c("CD14 Monocyte","CD16 Monocyte","Memory B","Naive B","Naive CD4 T","Naive CD8 T","Plasma","γδT")
my_comparison_order <- c("Mono_vs_Poly","Mono_vs_HC","Poly_vs_HC")

res <- plot_deg_pie_grid_summary(
  dat = summary_final2,
  celltype_order = my_celltype_order,
  comparison_order = my_comparison_order,
  col_up = "#EC706E",
  col_down = "#234091"
)

print(res$p_celltype_x)
print(res$p_comparison_x)




###(2) GSEA enrichment
library(dplyr)
library(clusterProfiler)
library(msigdbr)

run_deg_gsea <- function(de_res_filter, species = "Homo sapiens", pvalueCutoff = 0.2){
  
  genelist_db <- list(
    KEGG = msigdbr(species = species, category = "C2", subcategory = "CP:KEGG"),
    GO_BP = msigdbr(species = species, category = "C5", subcategory = "GO:BP"),
    REACTOME = msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME"),
    HALLMARK = msigdbr(species = species, category = "H")
  )
  genelist_db <- lapply(genelist_db, function(x) split(x$gene_symbol, x$gs_name))
  
  out_all <- list()
  cell_vec <- names(de_res_filter)
  
  for(ct in cell_vec){
    df_ct <- de_res_filter[[ct]]
    comp_vec <- unique(df_ct$comparison)
    
    for(comp in comp_vec){
      df_sub <- df_ct %>% filter(comparison == comp)
      if(nrow(df_sub) < 10) next
      
      gene_rank <- df_sub$avg_log2FC * (-log10(df_sub$p_val_adj + 1e-300))
      names(gene_rank) <- df_sub$gene
      gene_rank <- sort(gene_rank, decreasing = TRUE)
      
      for(db in names(genelist_db)){
        term2gene <- stack(genelist_db[[db]])[,c(2,1)]
        gsea_obj <- GSEA(geneList = gene_rank, TERM2GENE = term2gene, pvalueCutoff = pvalueCutoff, verbose = FALSE)
        df_gsea <- as.data.frame(gsea_obj@result)
        if(nrow(df_gsea) == 0) next
        
        df_gsea$celltype <- ct
        df_gsea$comparison <- comp
        df_gsea$database <- db
        df_gsea$Description <- gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_","",df_gsea$Description)
        out_all[[length(out_all)+1]] <- df_gsea
      }
    }
  }
  res_df <- bind_rows(out_all)
  return(res_df)
}

gsea_full_df <- run_deg_gsea(de_res_filter)
write.csv(gsea_full_df,"./CellPreprocess/gsea_all_celltype_comparison.csv",row.names=F)


###(3) Lolliplot 
library(readxl)
gsea_df <- read.csv("./CellPreprocess/gsea_all_celltype_comparison.csv",
                    stringsAsFactors = FALSE, check.names = FALSE)

View(gsea_df)

celltype_levels <- c("Naive CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                     "Tgd","Naive B","Memory B","Plasma","CD14 Monocyte",
                     "CD16 Monocyte","pDC","Platelet","Erythrocyte")
cell_color_vec <- c("#847AB3","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                    "#F1BAAB","#BB0021", "#D9C2D9","#FEC260","#AF478A",
                    "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")



library(tidyverse)
library(tidytext)

plot_gsea_lollipop <- function(gsea_df,
                               celltype_levels,
                               cell_color_vec,
                               size_range = c(1,5),
                               x_angle = 45){
  
  df_pre <- gsea_df %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    mutate(
      size_val = -log10(p.adjust + 1e-10),
      comparison = factor(comparison, levels = c("Poly_vs_HC","Mono_vs_Poly","Mono_vs_HC")),
      celltype = factor(celltype, levels = celltype_levels),
      x_id = str_c(celltype,"||",Description)
    )
  
  df_pre <- df_pre %>%
    group_by(comparison) %>%
    arrange(desc(NES > 0), desc(ifelse(NES > 0, NES, -Inf)), ifelse(NES < 0, abs(NES), Inf), .by_group = TRUE) %>%
    mutate(x_order = row_number()) %>%
    ungroup()
  
  label_df <- df_pre %>% distinct(x_id, Description)
  
  p <- ggplot(df_pre, aes(x = reorder_within(x_id, x_order, comparison), y = NES)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_segment(aes(xend = reorder_within(x_id, x_order, comparison), y = 0, yend = NES), colour = "gray50", linewidth = 0.3) +
    geom_point(aes(size = size_val, fill = celltype), shape = 21, stroke = 0) +
    scale_fill_manual(values = cell_color_vec) +
    scale_size_continuous(name = "-log10(adj.p)", range = size_range) +
    scale_x_reordered(labels = function(x){
      x_clean <- gsub("__.*$","",x)
      label_df$Description[match(x_clean, label_df$x_id)]
    }) +
    facet_wrap(~comparison, scales = "free_x", ncol = 3) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = x_angle, hjust = 1, size = 7, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_blank(),
      strip.text = element_text(colour = "black"),
      panel.spacing = unit(0.3,"cm")
    ) +
    labs(x = "Pathway", y = "NES")
  return(p)
}


p_gsealolli <- plot_gsea_lollipop(gsea_df, celltype_levels, 
        cell_color_vec, size_range = c(1,6.5), x_angle = 90)
p_gsealolli

###############################
###############################
################pseudobulk DEG
###### Load packages
library(Seurat); library(dplyr); library(tidyr); library(tibble)
library(DESeq2); library(clusterProfiler); library(msigdbr)
library(scatterpie); library(tidytext)

###### 1. Pseudobulk DEG
cell_types <- unique(pbmc1$celltype)
comp_list <- list(
  Poly_vs_HC = c("Polygenic JIA", "Healthy controls"),
  Mono_vs_Poly = c("Monogenic JIA", "Polygenic JIA"),
  Mono_vs_HC = c("Monogenic JIA", "Healthy controls")
)
deg_list <- list()
for(ct in cell_types){
  sub <- subset(pbmc1, celltype == ct)
  if(ncol(sub) < 10) next
  for(nm in names(comp_list)){
    g1 <- comp_list[[nm]][1]; g2 <- comp_list[[nm]][2]
    m <- sub@meta.data
    keep <- m$status %in% c(g1, g2)
    if(sum(keep) < 10) next
    sub2 <- sub[, keep]
    m2 <- sub2@meta.data
    if(!any(m2$status == g1) | !any(m2$status == g2)) next
    s1 <- unique(m2$sample[m2$status == g1])
    s2 <- unique(m2$sample[m2$status == g2])
    if(length(s1) < 2 | length(s2) < 2) next
    pseudo_mat <- AggregateExpression(sub2, group.by = "sample", assays = "RNA",
                                      slot = "counts", return.seurat = FALSE)$RNA
    if(is.null(pseudo_mat) | ncol(pseudo_mat) < 2) next
    pseudo_meta <- m2[!duplicated(m2$sample), c("sample","status")]
    rownames(pseudo_meta) <- pseudo_meta$sample
    common_samples <- intersect(colnames(pseudo_mat), rownames(pseudo_meta))
    if(length(common_samples) < 2) next
    pseudo_mat <- pseudo_mat[, common_samples, drop = FALSE]
    pseudo_meta <- pseudo_meta[common_samples, , drop = FALSE]
    n_g1 <- sum(pseudo_meta$status == g1); n_g2 <- sum(pseudo_meta$status == g2)
    if(n_g1 < 2 | n_g2 < 2) next
    pseudo_meta$sample <- NULL
    pseudo_meta$status <- factor(pseudo_meta$status, levels = c(g2, g1))
    if(nrow(pseudo_mat) == 0 | ncol(pseudo_mat) == 0) next
    if(sum(pseudo_mat) == 0) next
    pseudo_mat <- pseudo_mat[rowSums(pseudo_mat) > 0, , drop = FALSE]
    if(nrow(pseudo_mat) == 0) next
    dds <- DESeqDataSetFromMatrix(countData = round(pseudo_mat), colData = pseudo_meta, design = ~ status)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("status", g1, g2))
    res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>%
      mutate(comparison = nm, celltype = ct)
    deg_list[[paste0(ct, "_", nm)]] <- res_df
  }
}
deg_df <- bind_rows(deg_list) %>% filter(!is.na(padj))
saveRDS(deg_sig,file="E:/Cohort PPT/JIA/code/CellPreprocess/0920_degdf.RDS")
deg_sig <- deg_df %>% filter(abs(log2FoldChange) > 1, padj < 0.01)
deg_sig$direction <- ifelse(deg_sig$log2FoldChange > 0, "up", "down")
dim(deg_sig)
saveRDS(deg_sig,file="E:/Cohort PPT/JIA/code/CellPreprocess/0920_degsig.RDS")

deg_sig<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\0920_degsig.RDS")
deg_df<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\0920_degdf.RDS")


###### 2. Summary and gene set intersection
#####(1) summary
load("E:\\通路基因集合\\gene4pathway.Rdata")
IFN28 <- c("DDX60","EPSTI1","HERC5","HERC6","IFI27","IFI44","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFIT5","ISG15","LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SIGLEC1","SPATS2L","USP18","CXCL10","GBP1","SOCS1")
pyrop<-c("AIM2","APIP","CASP1","CASP4","CASP8","DHX9","ELANE","GSDMA","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB","NAIP","NLRC4","NLRP1","NLRP6","NLRP9","TREM2","ZBP1")#form Spectra
reactome <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:REACTOME")
IL1_genes <- reactome %>%filter(gs_name == "REACTOME_INTERLEUKIN_1_SIGNALING") %>%pull(gene_symbol) %>%unique()

gene_sets <- list(IFN28 = IFN28, IFN2 = IFN2$gene, MAPK = MAPK$gene,
                  NFKB = NFKB$gene, pyrop = pyrop, IL1 = IL1_genes)
stat_df <- deg_sig %>% group_by(celltype, comparison, direction) %>%
  summarise(gene_count = n(), .groups = "drop")
split_list <- split(deg_sig, paste(deg_sig$celltype, deg_sig$comparison, deg_sig$direction, sep = "_"))
intersect_res <- lapply(split_list, function(df){
  out <- lapply(gene_sets, function(s) length(intersect(df$gene, s)))
  tibble(!!!out, celltype = unique(df$celltype),
         comparison = unique(df$comparison), direction = unique(df$direction))
}) %>% bind_rows()
summary_final <- left_join(stat_df, intersect_res, by = c("celltype","comparison","direction"))


#####(2) DEG gene scatters
plot_manhattan_scatter <- function(deg_sig, celltypes_show,
                                   n_label = 15,
                                   bg_color = "gray90",
                                   bg_size = 0.5,
                                   pathway_size = 1.5,
                                   label_stroke = 0.3,
                                   text_size = 2,
                                   direction = "up",
                                   show_nfkb = TRUE, show_mapk = TRUE, show_ifn28 = TRUE,
                                   show_ifn2 = TRUE, show_pyrop = FALSE, show_il1 = FALSE,
                                   pathway_colors = c(
                                     NFKB = "#8C2522", MAPK = "#F49D5C",
                                     IFN28 = "#6699CC", IFN2 = "#B3D1E7",
                                     pyrop = "#A992C0", IL1 = "#AF478A"
                                   )){
  deg_vol <- deg_sig %>% filter(!is.na(padj), !is.na(log2FoldChange))
  all_genes <- list(
    NFKB = NFKB$gene, MAPK = MAPK$gene,
    IFN28 = IFN28, IFN2 = IFN2$gene,
    pyrop = pyrop, IL1 = IL1_genes
  )
  show_flags <- c(show_nfkb, show_mapk, show_ifn28, show_ifn2, show_pyrop, show_il1)
  keep_paths <- names(all_genes)[show_flags]
  col_use <- pathway_colors[keep_paths]
  man_list <- list()
  for(comp in c("Mono_vs_HC","Mono_vs_Poly","Poly_vs_HC")){
    df <- deg_vol %>%
      filter(comparison == comp, celltype %in% celltypes_show) %>%
      mutate(pathway = "background")
    for(pw in keep_paths){
      df$pathway[df$gene %in% all_genes[[pw]]] <- pw
    }
    if(direction == "up") df$pathway[df$pathway != "background" & df$log2FoldChange < 0] <- "background"
    if(direction == "down") df$pathway[df$pathway != "background" & df$log2FoldChange > 0] <- "background"
    df_sig <- df %>% filter(pathway != "background")
    top_genes <- df_sig %>%
      group_by(pathway) %>%
      arrange(desc(abs(log2FoldChange)), .by_group = TRUE) %>%
      slice_head(n = n_label) %>%
      pull(gene)
    df$label <- ifelse(df$gene %in% top_genes, df$gene, NA)
    df$is_labeled <- df$gene %in% top_genes
    df$celltype <- factor(df$celltype, levels = celltypes_show)
    df$pathway <- factor(df$pathway, levels = c("background", keep_paths))
    df_bg <- df %>% filter(pathway == "background")
    df_path_all <- df %>% filter(pathway != "background", !is_labeled)
    df_path_label <- df %>% filter(pathway != "background", is_labeled)
    p <- ggplot() +
      geom_jitter(data = df_bg, aes(x = celltype, y = log2FoldChange), width = 0.2, size = bg_size, color = bg_color, shape = 19) +
      geom_jitter(data = df_path_all, aes(x = celltype, y = log2FoldChange, color = pathway), width = 0.2, size = pathway_size, shape = 19) +
      geom_jitter(data = df_path_label, aes(x = celltype, y = log2FoldChange, fill = pathway), width = 0.2, size = pathway_size, shape = 21, stroke = label_stroke, color = "black") +
      geom_hline(yintercept = c(-1,1), linetype = "dashed", color = "grey50") +
      geom_text_repel(data = df_path_label, aes(x = celltype, y = log2FoldChange, label = label, color = pathway), size = text_size, max.overlaps = 30, na.rm = TRUE, show.legend = FALSE) +
      scale_color_manual(values = col_use) +
      scale_fill_manual(values = col_use) +
      coord_cartesian(ylim = c(-10, 10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
            axis.text.y = element_text(size = 7, color = "black"),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(x = "", y = "log2FoldChange", title = comp)
    man_list[[comp]] <- p
  }
  return(man_list)
}

celltypes_show <- c("CD14 Monocyte","CD16 Monocyte","pDC","Cytotoxic CD8 T",
                    "Naive CD8 T","Naive CD4 T","NK","Treg","γδT",
                    "Memory B","Naive B","Plasma")
man_list <- plot_manhattan_scatter(deg_sig, celltypes_show, direction = "up")
man_list$Mono_vs_HC
man_list$Mono_vs_Poly
man_list$Poly_vs_HC


#####(3) stacked plots
plot_pathway_stack <- function(summary_final, comp, celltypes_show,
                               pathway_order = c("NFKB","MAPK","IFN28","IFN2","pyrop","IL1"),
                               pathway_cols = c(NFKB = "#8C2522", MAPK = "#F49D5C",
                                                IFN28 = "#6699CC", IFN2 = "#B3D1E7",
                                                pyrop = "#A992C0", IL1 = "#AF478A")){
  df <- summary_final %>%
    filter(comparison == comp, celltype %in% celltypes_show) %>%
    select(celltype, direction, IFN28, IFN2, MAPK, NFKB, pyrop, IL1) %>%
    pivot_longer(cols = c(IFN28, IFN2, MAPK, NFKB, pyrop, IL1),
                 names_to = "pathway", values_to = "count") %>%
    filter(count > 0) %>%
    mutate(celltype = factor(celltype, levels = celltypes_show),
           direction = factor(direction, levels = c("up","down")),
           pathway = factor(pathway, levels = pathway_order))
  p <- ggplot(df, aes(x = direction, y = count, fill = pathway)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.2) +
    facet_wrap(~ celltype, nrow = 1) +
    scale_fill_manual(values = pathway_cols) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7, color = "black"),
          strip.text = element_text(size = 7, angle = 45),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "", y = "Gene count", title = comp, fill = "Pathway")
  return(p)
}

celltypes_show <- c("CD14 Monocyte","CD16 Monocyte","pDC","Cytotoxic CD8 T",
                    "Naive CD8 T","Naive CD4 T","NK","Treg","γδT",
                    "Memory B","Naive B","Plasma")
stack_p1 <- plot_pathway_stack(summary_final, "Mono_vs_HC", celltypes_show)
stack_p2 <- plot_pathway_stack(summary_final, "Poly_vs_HC", celltypes_show)
stack_p3 <- plot_pathway_stack(summary_final, "Mono_vs_Poly", celltypes_show)

sta_scatter1 <- (stack_p1+theme(legend.position = "none"))/ (man_list$Mono_vs_HC+theme(plot.title = element_blank())) + plot_layout(heights = c(0.5, 1.2))
sta_scatter2 <- (stack_p2+theme(legend.position = "none")) / (man_list$Poly_vs_HC+theme(plot.title = element_blank())) + plot_layout(heights = c(0.5, 1.2))
sta_scatter3 <- (stack_p3+theme(legend.position = "none"))/ (man_list$Mono_vs_Poly+theme(plot.title = element_blank())) + plot_layout(heights = c(0.5, 1.2))
ggsave(sta_scatter1, file = "E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\scatter_stack1.pdf", width = 7.95, height = 5.86)
ggsave(sta_scatter2, file = "E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\scatter_stack2.pdf", width = 7.95, height = 5.86)
ggsave(sta_scatter3, file = "E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\scatter_stack3.pdf", width = 7.95, height = 5.86)

sta_scatter1/sta_scatter2/sta_scatter3->sta
ggsave(sta, file = "E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\scatter_stack3.pdf", width = 7.95, height = 14.65)


#####(5) DEG count heatmap
plot_pathway_heatmap_grid <- function(deg_sig, celltypes_show,
                                      comparison_order = c("Mono_vs_HC","Mono_vs_Poly","Poly_vs_HC"),
                                      pathway_order = c("NFKB","MAPK","IFN28","IFN2","pyrop","IL1"),
                                      show_nfkb = TRUE, show_mapk = TRUE, show_ifn28 = TRUE,
                                      show_ifn2 = TRUE, show_pyrop = TRUE, show_il1 = TRUE,
                                      col_up = "#EC706E", col_down = "#234091",
                                      tile_w = 1, tile_h = 1){
  library(tidyverse)
  library(patchwork)
  deg_vol <- deg_sig %>% filter(!is.na(padj), !is.na(log2FoldChange))
  all_genes <- list(NFKB=NFKB$gene, MAPK=MAPK$gene, IFN28=IFN28,
                    IFN2=IFN2$gene, pyrop=pyrop, IL1=IL1_genes)
  show_flags <- c(show_nfkb, show_mapk, show_ifn28, show_ifn2, show_pyrop, show_il1)
  keep_paths <- pathway_order[show_flags]
  rows <- list()
  for(ct in celltypes_show){
    for(comp in comparison_order){
      sub <- deg_vol %>% filter(celltype == ct, comparison == comp)
      for(pw in keep_paths){
        g <- sub %>% filter(gene %in% all_genes[[pw]])
        n_up <- sum(g$log2FoldChange > 0, na.rm = TRUE)
        n_dn <- sum(g$log2FoldChange < 0, na.rm = TRUE)
        rows[[length(rows)+1]] <- data.frame(celltype=ct, comparison=comp, pathway=pw, direction="up", n=n_up)
        rows[[length(rows)+1]] <- data.frame(celltype=ct, comparison=comp, pathway=pw, direction="down", n=n_dn)
      }
    }
  }
  df <- do.call(rbind, rows)
  up_max <- max(df$n[df$direction=="up"])
  dn_max <- max(df$n[df$direction=="down"])
  df <- df %>% mutate(
    val = ifelse(direction=="up", n/up_max, -n/dn_max),
    xlab = paste0(pathway, ifelse(direction=="up","\u2191","\u2193"))
  )
  x_levels <- as.vector(sapply(keep_paths, function(pw) paste0(pw, c("\u2191","\u2193"))))
  df$xlab <- factor(df$xlab, levels=x_levels)
  plots <- list()
  for(i in seq_along(celltypes_show)){
    ct <- celltypes_show[i]
    sub <- df %>% filter(celltype == ct)
    p <- ggplot(sub, aes(x=xlab, y=comparison, fill=val)) +
      geom_tile(color="white", linewidth=0.6, width=tile_w, height=tile_h) +
      geom_text(aes(label=n), size=2.4, color="black") +
      scale_fill_gradient2(low=col_down, mid="white", high=col_up,
                           midpoint=0, limits=c(-1,1), breaks=c(-1,0,1),
                           labels=c(paste0(up_max,"\u2193"),"0",paste0(up_max,"\u2191")), name="Gene count") +
      scale_y_discrete(limits=comparison_order) +
      labs(y=ct) +
      theme_classic() +
      theme(axis.text.x = element_text(size=7, color="black"),#angle=45, hjust=1, 
            axis.text.y = element_text(size=7, color="black"),
            axis.title.y = element_text(size=8, color="black"),
            axis.title.x = element_blank(),
            legend.position="none",
            plot.margin = margin(0,0,0,0,"pt"))
    if(i < length(celltypes_show)){
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    plots[[i]] <- p
  }
  combined <- wrap_plots(plots, ncol=1)
  return(combined)
}

celltypes_show <- c("CD14 Monocyte","CD16 Monocyte","pDC","Cytotoxic CD8 T",
                    "Naive CD8 T","Naive CD4 T","NK","Treg","γδT",
                    "Memory B","Naive B","Plasma")
p_heat <- plot_pathway_heatmap_grid(deg_sig, celltypes_show,
                                    col_up = "#D48A47", col_down = "#234091")
p_heat



###### 3. GSEA, single cell type, single sample
library(msigdb)
library(msigdbr)
library(clusterProfiler)
library(Seurat)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(msigdbr)

run_gsea_single <- function(obj, celltypes, mutant_samples, hc_samples, poly_samples,
                            species = "Homo sapiens", pvalueCutoff = 0.2){
  options(future.globals.maxSize = Inf)
  future::plan(future::sequential)
  hallmark <- msigdbr(species = species, category = "H")
  term2gene <- stack(split(hallmark$gene_symbol, hallmark$gs_name))[, c(2,1)]
  out <- list()
  for(ct in celltypes){
    sub_ct <- subset(obj, subset = celltype == ct)
    for(mut in mutant_samples){
      for(ref_grp in list(list(ref = hc_samples, tag = "vs_HC"),
                          list(ref = poly_samples, tag = "vs_Poly"))){
        options(future.globals.maxSize = Inf)
        future::plan(future::sequential)
        keep_samples <- c(mut, ref_grp$ref)
        sub <- subset(sub_ct, subset = sample %in% keep_samples)
        sub$group <- ifelse(sub$sample == mut, "mut", "ref")
        Idents(sub) <- sub$group
        markers <- FindMarkers(sub, ident.1 = "mut", ident.2 = "ref",
                               logfc.threshold = 0, min.pct = 0, verbose = FALSE)
        markers <- markers %>% tibble::rownames_to_column("gene")
        if(nrow(markers) < 10) next
        markers$score <- markers$avg_log2FC * (-log10(markers$p_val_adj + 1e-300))
        gene_rank <- sort(markers$score, decreasing = TRUE)
        names(gene_rank) <- markers$gene
        gsea_obj <- GSEA(geneList = gene_rank, TERM2GENE = term2gene,
                         pvalueCutoff = pvalueCutoff, verbose = FALSE)
        if(nrow(gsea_obj@result) == 0) next
        df_gsea <- as.data.frame(gsea_obj@result)
        df_gsea$Description <- gsub("HALLMARK_", "", df_gsea$Description)
        df_gsea$celltype <- ct
        df_gsea$mutant <- mut
        df_gsea$comparison <- paste0(mut, ref_grp$tag)
        out[[length(out)+1]] <- df_gsea
      }
    }
  }
  res_df <- bind_rows(out)
  return(res_df)
}

mutant_samples <- c("182LPIN2","181NOD2","183NOD2","190PSTPIP1","200PSTPIP1")
hc_samples <- paste0("C", 1:6)
poly_samples <- paste0("poly", 1:6)
celltypes_use <- c("CD14 Monocyte","CD16 Monocyte","pDC","Cytotoxic CD8 T")
celltypes_use2<-c("NK","Plasma","Naive CD4 T","Naive CD8 T")

gsea_res <- run_gsea_single(pbmc1, celltypes_use, mutant_samples, hc_samples, poly_samples)
gsea_res2 <- run_gsea_single(pbmc1, celltypes_use, mutant_samples, hc_samples, poly_samples)


dim(gsea_res)
head(gsea_res)

openxlsx::write.xlsx(gsea_res,file="E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\0923_GSEA.xlsx")
openxlsx::write.xlsx(gsea_res2,file="E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\0923_GSEA2.xlsx")

###bubble plot (p.adjust<0.05)
multiplesheets <- function(fname) {
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  names(data_frame) <- sheets
  print(data_frame)
}

### monogenic VS HC
plot_gsea_bubble <- function(data, fill_palette, pathway_order, pathway_labels, shape_type = 16, point_scale = 1.2){
  library(ggplot2); library(dplyr)
  data <- data[data$type == "VSHC", ]
  data$NES <- as.numeric(data$NES)
  data$padj <- as.numeric(data$p.adjust)
  data$neglog10padj <- -log10(data$padj)
  data$ID <- factor(data$ID, levels = pathway_order, labels = pathway_labels)
  mutant_order <- c("182LPIN2", "181NOD2", "183NOD2", "200PSTPIP1", "190PSTPIP1")
  mutant_labels <- c("182LPIN2" = "P9", "181NOD2" = "P13", "183NOD2" = "P14",
                     "200PSTPIP1" = "P15", "190PSTPIP1" = "P16")
  data$mutant <- factor(data$mutant, levels = mutant_order, labels = mutant_labels[mutant_order])
  data <- droplevels(data)
  base_theme <- theme_bw() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          strip.text = element_text(color = "black", size = 9),
          legend.title = element_text(color = "black", size = 8),
          legend.text = element_text(color = "black", size = 7),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.5))
  p_col <- ggplot(data, aes(x = ID, y = celltype, size = NES, color = neglog10padj)) +
    geom_point(shape = shape_type, stroke = 0.5) +
    scale_color_gradientn(colors = fill_palette) +
    scale_size_continuous(range = c(2*point_scale, 8*point_scale)) +
    facet_grid(mutant ~ ., scales = "free_y", space = "free_y") +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "", size = "NES", color = "-log10(p.adjust)")
  p_row <- ggplot(data, aes(x = celltype, y = ID, size = NES, color = neglog10padj)) +
    geom_point(shape = shape_type, stroke = 0.5) +
    scale_color_gradientn(colors = fill_palette) +
    scale_size_continuous(range = c(2*point_scale, 8*point_scale)) +
    facet_grid(. ~ mutant, scales = "free_x", space = "free_x") +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "", size = "NES", color = "-log10(p.adjust)")
  return(list(column = p_col, row = p_row))
}
pathway_order <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING")
pathway_labels <- c("Interferon-γ response", "Interferon-α response",
                    "Inflammatory response", "TNF via NF-κB", "IL6-JAK-STAT3")
col2<-colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
resu<-multiplesheets("E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\0923_GSEA.xlsx")
resu$`VS HC and poly`->data
data[data$type=="VSHC",]->data1
p <- plot_gsea_bubble(data1, col2, pathway_order, pathway_labels, shape_type = 15, point_scale = 1.2)
print(p)


### monogenic VS polygenic
data2 <- resu$Polybar
plot_gsea_bubble <- function(data, fill_palette, pathway_order, pathway_labels,
                             shape_type = 16, point_scale = 1.2){
  library(ggplot2); library(dplyr)
  data$NES <- as.numeric(data$NES)
  data$padj <- as.numeric(data$p.adjust)
  data$neglog10padj <- -log10(data$padj)
  data$ID <- factor(data$ID, levels = pathway_order, labels = pathway_labels)
  mutant_order <- c("182LPIN2", "181NOD2", "183NOD2", "200PSTPIP1", "190PSTPIP1")
  mutant_labels <- c("182LPIN2" = "P9", "181NOD2" = "P13", "183NOD2" = "P14",
                     "200PSTPIP1" = "P15", "190PSTPIP1" = "P16")
  data$mutant <- factor(data$mutant, levels = mutant_order, labels = mutant_labels[mutant_order])
  data <- droplevels(data)
  base_theme <- theme_bw() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          strip.text = element_text(color = "black", size = 9),
          legend.title = element_text(color = "black", size = 8),
          legend.text = element_text(color = "black", size = 7),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.5))
  p_col <- ggplot(data, aes(x = ID, y = celltype, size = NES, color = neglog10padj)) +
    geom_point(shape = shape_type, stroke = 0.5) +
    scale_color_gradientn(colors = fill_palette) +
    scale_size_continuous(range = c(2*point_scale, 8*point_scale)) +
    facet_grid(mutant ~ ., scales = "free_y", space = "free_y") +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "", size = "NES", color = "-log10(p.adjust)")
  p_row <- ggplot(data, aes(x = celltype, y = ID, size = NES, color = neglog10padj)) +
    geom_point(shape = shape_type, stroke = 0.5) +
    scale_color_gradientn(colors = fill_palette) +
    scale_size_continuous(range = c(2*point_scale, 8*point_scale)) +
    facet_grid(. ~ mutant, scales = "free_x", space = "free_x") +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "", size = "NES", color = "-log10(p.adjust)")
  return(list(column = p_col, row = p_row))
}

pathway_order <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING")
pathway_labels <- c("Interferon-γ response", "Interferon-α response",
                    "Inflammatory response", "TNF via NF-κB", "IL6-JAK-STAT3")
col2 <- colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
p_poly <- plot_gsea_bubble(data2, col2, pathway_order, pathway_labels,
                           shape_type = 15, point_scale = 1.2)
p_poly$column
p_poly$row
