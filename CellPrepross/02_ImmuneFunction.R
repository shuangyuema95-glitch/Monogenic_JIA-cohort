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
  stats_list <- list()
  for(ct in cell_order){
    ct_data <- DATA[DATA$celltype == ct, ]
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
  if(length(stats_list) == 0){
    stats_df <- data.frame(celltype=character(), group1=character(), group2=character(),
                           n1=numeric(), n2=numeric(), median1=numeric(), median2=numeric(),
                           p_value=numeric(), p_adj=numeric(), stringsAsFactors = FALSE)
    sig_df <- stats_df
  } else {
    stats_df <- do.call(rbind, stats_list) %>% group_by(celltype) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>% arrange(celltype, p_adj) %>% as.data.frame()
    sig_df <- stats_df[stats_df$p_adj < 0.05, ]
  }
  plot_list <- list()
  for(i in seq_along(cell_order)){
    ct <- cell_order[i]
    ct_data <- DATA[DATA$celltype == ct, c("group", score)]
    colnames(ct_data)[2] <- "score"
    ct_data <- ct_data[!is.na(ct_data$score), ]
    hc_mean <- mean(ct_data$score[ct_data$group == "HC"], na.rm = TRUE)
    p <- ggplot(ct_data, aes(x = group, y = score, fill = group)) +
      geom_violin(trim = FALSE, scale = "width", width = violin_width, color = NA) +
      geom_boxplot(width = boxplot_width, color = "black", outlier.shape = NA, fill = NA) +
      geom_hline(yintercept = hc_mean, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      scale_fill_manual(values = all_colors) +
      labs(title = ct, x = "", y = ifelse(i == 1, score, "")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            plot.margin = margin(2, 2, 2, 2), legend.position = "none",
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.y = element_line(color = "black"),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.y = element_text(colour = "black", size = 9))
    plot_list[[ct]] <- p
  }
  combined <- wrap_plots(plot_list, ncol = length(cell_order)) + plot_layout(guides = "collect") &
    theme(plot.margin = margin(2, 2, 2, 2))
  cat("Score:", score, "| Cells:", paste(cell_order, collapse = ", "),
      "| Total comparisons:", nrow(stats_df),
      "| Significant (adj<0.05):", nrow(sig_df), "\n")
  return(list(plot = combined, stats = stats_df, sig = sig_df))
}

res_nfkb <- Violin_All_row(pbmc1, score = "auc_NFKB",violin_width = 0.65, boxplot_width = 0.68)
res_nfkb$plot
res_mapk <- Violin_All_row(pbmc1, score = "auc_MAPK",violin_width = 0.65, boxplot_width = 0.68)
res_mapk$plot
res_ifn2 <- Violin_All_row(pbmc1, score = "auc_IFNII",violin_width = 0.65, boxplot_width = 0.68)
res_ifn2$plot



all_cells <- unique(pbmc1$celltype)
res_nfkb <- Violin_All_row(pbmc1, score = "auc_NFKB", cell_order = all_cells,violin_width = 0.65, boxplot_width = 0.68)
res_nfkb$plot
res_mapk <- Violin_All_row(pbmc1, score = "auc_MAPK", cell_order = all_cells,violin_width = 0.65, boxplot_width = 0.68)
res_mapk$plot
res_ifn2 <- Violin_All_row(pbmc1, score = "auc_IFNII", cell_order = all_cells,violin_width = 0.65, boxplot_width = 0.68)
res_ifn2$plot



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



