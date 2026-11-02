library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(ggplot2)
library(DoubletFinder)
library(patchwork)
library(cowplot)
library(openxlsx)

setwd("E:/Cohort PPT/JIA/code")

samples=c("C1","C2","C3","C4","C5","C6","poly1","poly2","poly3","poly4","poly5","poly6","157LACC1","182LPIN2","181NOD2","183NOD2","190PSTPIP1","200PSTPIP1")
gene_variants<-c(rep("none gene",12),c("LACC1_S366*/N19Kfs*31","LPIN2_I67R/K781Tfs*7","NOD2_M513K","NOD2_C495Y","PSTPIP1_E250K","PSTPIP1_E250K"))
length(samples)==length(gene_variants)
fold2="E:/Cohort PPT/JIA/code"


#####1 ----preprocessing-----
sce <- list()
summary_df <- data.frame()
for(i in 1:length(samples)) {
  s <- CreateSeuratObject(counts = Read10X(paste(fold2, samples[i], sep="/")),
                          min.cells = 3, min.features = 200)
  
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s[["percent.ribo"]] <- PercentageFeatureSet(s, pattern = "^RP[SL]")
  s$log10GenesPerUMI <- log10(s$nFeature_RNA) / log10(s$nCount_RNA)
  
  s <- subset(s,
              subset = nCount_RNA > 1000 & nCount_RNA < 100000 &
                nFeature_RNA > 500 & nFeature_RNA < 8000 &
                percent.mt < 20 &
                log10GenesPerUMI > 0.7)
  n_qc <- ncol(s)
  
  all.genes <- rownames(s)
  exclude.genes <- unique(c(
    grep("^RPS", all.genes, value = TRUE), # we excluded ribosome genes
    grep("^RPL", all.genes, value = TRUE), # we excluded ribosome genes
    grep("^HSP", all.genes, value = TRUE), # we excluded heat shock genes
    grep("^MT-", all.genes, value = TRUE)  # we excluded mitchondrial gens
  ))
  keep.genes <- setdiff(all.genes, exclude.genes)
  
  # Normalization and FindVariable for each datasets and before removing doublets
  s <- subset(s, features = keep.genes)
  s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000)
  s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 4000)
  hvg_genes <- VariableFeatures(s)
  
  
  # Doublet testing
  npc <- min(50, ncol(s) - 1)
  pK <- paramSweep(s, PCs = 1:npc, sct = FALSE)
  sweep.stats <- summarizeSweep(pK, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  if (length(best.pK) == 0 || is.na(best.pK)) best.pK <- 0.01
  nExp <- round(0.03 * ncol(s))
  
  s <- doubletFinder(s, PCs = 1:npc, pN = 0.25, pK = best.pK,
                     nExp = nExp, reuse.pANN = FALSE, sct = FALSE)
  colnames(s@meta.data)[8] <- "scDbl_class"
  s <- subset(s, scDbl_class == "Singlet")
  
  
  # ---- annotate & collect ----
  s$sample <- samples[i]
  s$mut <- gene_variants[i]
  sce[[i]] <- s
  summary_df <- rbind(summary_df, c(samples[i], n_qc, ncol(s)))
}

names(sce)<-samples
summary_df
write.xlsx(summary_df,file="00_scCellcountDetails.xlsx")
saveRDS(sce,file="2_readsceList.RDS")

#####2 ----merge samples-----
sce<-readRDS("2_readsceList.RDS")
pbmc1 <- merge(x = sce[[1]], y = sce[-1], add.cell.ids = names(sce), project = "pbmc1")
table(pbmc1$sample)
table(sapply(rownames(pbmc1@meta.data), function(x){strsplit(x, split = "_")[[1]][1]}))
table(pbmc1$scDbl_class) #ncol is 174,561
#pbmc1 <- subset(pbmc1, scDbl_class == "Singlet")
###whole data normalization
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)

###check hvgs list intersect
# hvg_list <- lapply(sce, VariableFeatures)
# length(hvg_list )
# lapply(hvg_list ,length)
# hvgs <- Reduce("intersect",hvg_list)
# length(hvgs) #317
# hvgs <- hvgs[hvgs %in% rownames(pbmc1)]
# length(hvgs)
###check hvgs union
hvg_list <- lapply(sce, VariableFeatures)
hvgs <- Reduce("union",hvg_list)
hvgs <- hvgs[hvgs %in% rownames(pbmc1)]
length(hvgs)
gene_count <- sapply(hvgs, function(g){sum(sapply(hvg_list, function(hvgs_sample) g %in% hvgs_sample))})
threshold <- sort(gene_count, decreasing = TRUE)[1500]
hvgs_top1500 <- names(gene_count[gene_count >= threshold])
length(hvgs_top1500)#1789
hvgs<-hvgs_top1500

###correct 
VariableFeatures(pbmc1) <- hvgs
length(pbmc1@assays$RNA@var.features)
g2m_genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(pbmc1))
s_genes  <- CaseMatch(search = cc.genes$s.genes,  match = rownames(pbmc1))
pbmc1 <- CellCycleScoring(pbmc1, g2m.features = g2m_genes, s.features = s_genes)
pbmc1 <- ScaleData(pbmc1, vars.to.regress = c("percent.mt"))#, "S.Score", "G2M.Score"))
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(pbmc1), npcs = 50)
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(datasets = str_split_fixed(rownames(pbmc1@meta.data), "_", 2)[,1])
pbmc1 <- RunHarmony(pbmc1, group.by.vars = "datasets", dims.use = 1:50)# theta = 2) 


#####3 ----Cell clustering-----
pbmc1 <- FindNeighbors(pbmc1, dims = 1:50, reduction = "harmony")
pbmc1 <- FindClusters(pbmc1, resolution = 0.25)
length(unique(pbmc1$seurat_clusters))
set.seed(123)
pbmc1 <- RunUMAP(pbmc1, reduction = "harmony", dims = 1:50, 
                 reduction.name = "harmony_umap")
                # min.dist = 0.05, # default 0.00
                #n.neighbors = 50 #default 30)
pbmc1$harmony_clusters <- pbmc1$seurat_clusters

p1<-DimPlot(pbmc1, label = TRUE, reduction = "harmony_umap", group.by = "datasets",
            pt.size = 0.0001, raster = FALSE) + theme(plot.title = element_blank())
p2<-DimPlot(pbmc1, label = TRUE, reduction = "harmony_umap", split.by = "datasets",
            pt.size = 0.0001, raster = FALSE) + theme(plot.title = element_blank())
p3<-DimPlot(pbmc1, label = TRUE, reduction = "harmony_umap", group.by = "harmony_clusters",
            pt.size = 0.0001, raster = FALSE) + theme(plot.title = element_blank())

p1
p2
p3

#saveRDS(pbmc1,file="pbmc1.RDS")
pbmc1<-readRDS("pbmc1.RDS")



#####4 ----Cell annotation-----
library(COSG)
marker_cosg <- cosg(
  pbmc1,groups='all',assay='RNA',
  slot='data',mu=10,n_genes_user=100,
  remove_lowly_expressed=TRUE,
  expressed_pct=0.1)

View(marker_cosg$names)

pbmc1<-subset(pbmc1,seurat_clusters%in%c(1,2,3,4,5,6,7,8,9,10,
                                         11,12,13,14,15,16,17,18,20,
                                         22,23,24,25))

Idents(pbmc1) <- "seurat_clusters"
pbmc1<-RenameIdents(pbmc1,
                    '0'='Naive CD4 T',#"CD4","SELL","CCR7","TCF7" done 
                    '1'='Naive CD8 T',#"CD8A","CD8B","CCR7","LEF1"  done
                    '2'='CD14 Monocyte',#"CD14","S100A8","S100A9" done
                    '3'='Naive B',#"IGHD","IGHD","TCL1A","IL4R" done
                    '4'='NK', #"NKG7","KLRF1","KLRD1" done         
                    '5'='Cytotoxic CD8 T',# "CD8A","CD8B","GZMA","GZMB","GZMH" done 
                    '6'='Memory CD4 T',#"AQP3","TSHZ2","CD40LG","IL7R" done
                    '7' = 'γδT',#"TRGC1", "TRGC2","CXCR6"  done
                    '8' = 'Treg',#"FOXP3","IL2RA","IKZF2" done
                    '9'='Memory B',#"FCRL5", "FCRL2", "CD24","TNFRSF13B","CD80","IGHG3" done 
                    '10'='CD14 Monocyte',# #"CD14","S100A8","S100A9" done
                    '11'='Platelet',# "PF4","PPBP","ITGA2B" done
                    '12'='γδT',#"TRGC1", "TRGC2","CXCR6"  done
                    '13'='CD16 Monocyte',#"FCGR3A","MS4A7","LST1"  done        
                    '14'='Erythrocyte',#"HBB","HBA2","HBA1" done 
                    '15'='Naive B',#"IGHD","IGHD","TCL1A","IL4R" done 
                    '16'='NK', #"NKG7","KLRF1","KLRD1" done 
                    '17'='pDC',#"PLD4","IL3RA","IRF7","IL3RA" done 
                    '18'='CD14 Monocyte',# #"CD14","S100A8","S100A9" done
                    '19'='Plasma',#"JCHAIN","XBP1","MZB1", done CD19- CD38+ CD27+ 
                    '20'='Naive B',#"IGHD","IGHD","TCL1A","IL4R" done
                    #'21'='Prolif cell',#"TOP2A","MKI67"
                    '22'='CD14 Monocyte',# #"CD14","S100A8","S100A9" done
                    '23'='Naive B',#"IGHD","IGHD","TCL1A","IL4R" done
                    '24'='Cytotoxic CD8 T',# "CD8A","CD8B","GZMA","GZMB","GZMH" done 
                    '25'='Naive CD8 T'#"CD8A","CD8B","CCR7","LEF1"  done
)
pbmc1@meta.data$celltype <- Idents(pbmc1)

length(unique(pbmc1$celltype))#15
sample_map <- tibble(
  sample = c("C1","C2","C3","C4","C5","C6",
             "poly1","poly2","poly3","poly4","poly5","poly6",
             "157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1"),
  sample2 = c("C1","C2","C3","C4","C5","C6",
              "PJ1","PJ2","PJ3","PJ4","PJ5","PJ6",
             "P7", "P10","P14","P15","P17","P18"))
pbmc1@meta.data <- pbmc1@meta.data %>%left_join(sample_map, by = "sample")

pbmc1@meta.data<-pbmc1@meta.data%>%mutate(status=case_when(
  sample%in%c("C1","C2","C3","C4","C5","C6")~"Healthy controls",
  sample%in%c("poly1","poly2","poly3","poly4","poly5","poly6")~"Polygenic JIA",
  sample%in%c("157LACC1","182LPIN2","181NOD2","183NOD2","190PSTPIP1","200PSTPIP1")~"Monogenic JIA"))
status_level<-c("Healthy controls","Polygenic JIA","Monogenic JIA")

#####5 ----Cell UMAP-----
celltype_levels <- c("Naive CD4 T","Memory CD4 T",
                     "Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                     "γδT","Naive B","Memory B","Plasma","CD14 Monocyte",
                     "CD16 Monocyte","pDC","Platelet","Erythrocyte")
cell_color_vec <- c("#847AB3","#96AF95","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                    "#F1BAAB","#BB0021", "#D9C2D9","#FEC260","#AF478A",
                    "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")
names(cell_color_vec) <- celltype_levels

sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
cell_color_vec2 <- c("#847AB3","#96AF95","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                     "#F1BAAB","#BB0021", "#D9C2D9","#FEC260","#AF478A",
                     "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9","#B1DA99",
                     "#7B1FA2","#FFEFC1")
names(cell_color_vec2)<-sample_level

cell_color_vec3 <- c("#6699CC","#FEC260","#AF478A")
names(cell_color_vec3)<-status_level
table(pbmc1$sample)
table(pbmc1$status)
plot_umap_standard <- function(Seurat_obj,my_colors,reduction = "",group_by = "samples",pt_size = 0.15,pname="",leg_title=NULL,label_cluster = TRUE) {
  library(ggplot2)
  library(ggrepel)
  df <- as.data.frame(Seurat_obj@reductions[[reduction]]@cell.embeddings)
  colnames(df) <- c("UMAP_1","UMAP_2")
  df$Cluster <- as.factor(Seurat_obj@meta.data[[group_by]])
  use_title <- if(!is.null(leg_title)) leg_title else group_by
  centroid_df <- NULL
  if(label_cluster){
    centroid_df <- df %>% group_by(Cluster) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
  }
  p <- ggplot(df, aes(x=UMAP_1,y=UMAP_2,color=Cluster)) +
    geom_point(size=pt_size,shape=16,stroke=0) +
    scale_color_manual(values=my_colors) +
    theme_classic() +
    labs(color=use_title) +
    theme(plot.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.3),
          aspect.ratio=1,plot.title=element_blank(),
          legend.text=element_text(size=6),legend.title=element_text(size=6)) +
    guides(color=guide_legend(ncol=1,override.aes=list(size=2.5)))
  if(label_cluster){
    p <- p + theme(legend.position = "none") + geom_text_repel(data = centroid_df,aes(label = Cluster, color = Cluster),size = 3, max.overlaps = Inf,show.legend = FALSE)
  }else{
    p <- p + theme(legend.position="right")
  }
  return(p)
}

P1_celltype <- plot_umap_standard(pbmc1, 
                                  my_colors=cell_color_vec,reduction = "harmony_umap",
                                  group_by = "celltype", pt_size = 0.1, label_cluster = FALSE)

P1_samples <- plot_umap_standard(pbmc1, 
                                 my_colors=cell_color_vec2,reduction = "harmony_umap",
                                 group_by = "sample", pt_size = 0.1, label_cluster = FALSE)

P1_status <- plot_umap_standard(pbmc1, 
                                my_colors=cell_color_vec3,reduction = "harmony_umap",
                                group_by = "status", pt_size = 0.1, label_cluster = FALSE)

library(patchwork)
p1<-P1_celltype|P1_samples|P1_status
ggsave(p1,file="E:/Cohort PPT/JIA/code/00_0.05umap.pdf",width = 19.97,height = 3.25)
saveRDS(pbmc1,file="pbmc1.RDS")


#####6 ----Cell Proportion-----
####(1) total cell count to calculated cell fraction
PlotCountDistribution <- function(seu, Count_var, color_vec){
  library(purrr)
  df <- table(seu@meta.data[[Count_var]]) %>%
    as.data.frame() %>%
    set_names(c("group","n")) %>%
    mutate(frac = n / sum(n) * 100, label = sprintf("%.2f%%", frac)) %>%
    arrange(desc(frac))
  df$group <- factor(df$group, levels = rev(df$group))
  
  p_bar_h <- ggplot(df, aes(y = group, x = frac, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label), hjust = -0.15, size = 3, colour = "black") +
    scale_fill_manual(values = color_vec) +
    labs(x = "Percentage (%)", y = NULL) +
    theme_classic() +
    theme(legend.position = "none",axis.text = element_text(colour = "black"),plot.margin = margin(5,15,5,5)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, max(df$frac)*1.12))
  
  p_lollipop <- ggplot(df, aes(y = group, x = frac, color = group)) +
    geom_segment(aes(x = 0, xend = frac, y = group, yend = group), color = "grey60", linewidth = 0.8) +
    geom_point(size = 2.45) +
    geom_text(aes(label = label), hjust = -0.15, size = 3, colour = "black") +
    scale_color_manual(values = color_vec) +
    labs(x = "Percentage (%)", y = NULL) +
    theme_classic() +
    theme(legend.position = "none",axis.text = element_text(colour = "black"),plot.margin = margin(5,15,5,5)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, max(df$frac)*1.12))
  
  p_donut <- ggplot(df, aes(x = 2, y = frac, fill = group)) +
    geom_col(width = 1) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3, colour = "black") +
    scale_fill_manual(values = color_vec) +
    coord_polar("y", start = 0) +
    xlim(0.5,2.5) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(legend.title = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank())
  
  return(list(bar_h = p_bar_h, lollipop = p_lollipop, donut = p_donut))
}
plist_celltype <- PlotCountDistribution(pbmc1, Count_var = "celltype", color_vec = cell_color_vec)
plist_status <- PlotCountDistribution(pbmc1, Count_var = "status", color_vec = cell_color_vec3)
plist_celltype$lollipop/plist_status$donut
plist_celltype$donut/plist_status$donut
p_comb <- plist_celltype$lollipop + inset_element(
plist_status$donut+theme(legend.position = "none"),
left = 0.32, bottom = 0.018, right = 0.98, top = 0.38);p_comb 

PlotCountDistribution <- function(seu, Count_var, color_vec, split_by = NULL, split_keep = NULL, group_order = NULL){
  meta <- seu@meta.data
  if (!is.null(split_by)) {
    grp <- as.character(meta[[Count_var]])
    sub <- as.character(meta[[split_by]])
    grp[grp %in% split_keep] <- sub[grp %in% split_keep]
    meta[[Count_var]] <- grp
  }
  df <- table(meta[[Count_var]]) %>%
    as.data.frame() %>%
    set_names(c("group","n")) %>%
    mutate(frac = n / sum(n) * 100, label = sprintf("%.2f%%", frac)) %>%
    arrange(desc(frac))
  if (is.null(group_order)) df$group <- factor(df$group, levels = rev(df$group))
  else df$group <- factor(df$group, levels = group_order)
  p_bar_h <- ggplot(df, aes(y = group, x = frac, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label), hjust = -0.15, size = 3, colour = "black") +
    scale_fill_manual(values = color_vec) +
    labs(x = "Percentage (%)", y = NULL) +
    theme_classic() +
    theme(legend.position = "none", axis.text = element_text(colour = "black"), plot.margin = margin(5, 15, 5, 5)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(df$frac) * 1.12))
  p_lollipop <- ggplot(df, aes(y = group, x = frac, color = group)) +
    geom_segment(aes(x = 0, xend = frac, y = group, yend = group), color = "grey60", linewidth = 0.8) +
    geom_point(size = 2.45) +
    geom_text(aes(label = label), hjust = -0.15, size = 3, colour = "black") +
    scale_color_manual(values = color_vec) +
    labs(x = "Percentage (%)", y = NULL) +
    theme_classic() +
    theme(legend.position = "none", axis.text = element_text(colour = "black"), plot.margin = margin(5, 15, 5, 5)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(df$frac) * 1.12))
  p_donut <- ggplot(df, aes(x = 2, y = frac, fill = group)) +
    geom_col(width = 1) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3, colour = "black") +
    scale_fill_manual(values = color_vec) +
    coord_polar("y", start = 0) +
    xlim(0.5, 2.5) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(legend.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
  return(list(bar_h = p_bar_h, lollipop = p_lollipop, donut = p_donut))
}

PlotNestedDonut <- function(seu, outer_var, inner_var, outer_color, inner_color,
                            outer_level = NULL, inner_level = NULL,
                            outer_label = TRUE, inner_label = TRUE){
  mk <- function(v, lv) {
    d <- table(seu@meta.data[[v]]) %>% as.data.frame() %>% set_names(c("group","n")) %>%
      mutate(frac = n / sum(n) * 100, label = sprintf("%.2f%%", frac))
    if (!is.null(lv)) d$group <- factor(d$group, levels = lv)
    d
  }
  do <- mk(outer_var, outer_level)
  di <- mk(inner_var, inner_level)
  p <- ggplot() +
    geom_col(data = do, aes(x = 2.5, y = frac, fill = group), width = 1.4, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = outer_color) +
    ggnewscale::new_scale_fill() +
    geom_col(data = di, aes(x = 1.1, y = frac, fill = group), width = 1.0, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = inner_color) +
    coord_polar("y", start = 0) +
    xlim(0, 3.4) +
    theme_void() +
    theme(legend.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
  if (outer_label)
    p <- p + geom_text(data = do, aes(x = 2.5, y = frac, label = label),
                       position = position_stack(vjust = 0.5), size = 3, colour = "black")
  if (inner_label)
    p <- p + geom_text(data = di, aes(x = 1.1, y = frac, label = label),
                       position = position_stack(vjust = 0.5), size = 2.6, colour = "black")
  p
}

pbmc1@meta.data <- pbmc1@meta.data %>% mutate(genome = case_when(
  sample %in% c("C1","C2","C3","C4","C5","C6") ~ "Healthy controls",
  sample %in% c("poly1","poly2","poly3","poly4","poly5","poly6") ~ "Polygenic JIA",
  sample %in% c("157LACC1") ~ "LACC1",
  sample %in% c("182LPIN2") ~ "LPIN2",
  sample %in% c("181NOD2","183NOD2") ~ "NOD2",
  sample %in% c("190PSTPIP1","200PSTPIP1") ~ "PSTPIP1"))
genome_level <- c("Healthy controls","Polygenic JIA","NOD2","PSTPIP1","LPIN2","LACC1")
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(pathway = case_when(
  sample %in% c("C1","C2","C3","C4","C5","C6") ~ "Healthy controls",
  sample %in% c("poly1","poly2","poly3","poly4","poly5","poly6") ~ "Polygenic JIA",
  sample %in% c("157LACC1") ~ "Immune metabolism",
  sample %in% c("182LPIN2") ~ "Inflammasome",
  sample %in% c("181NOD2","183NOD2") ~ "NFKB",
  sample %in% c("190PSTPIP1","200PSTPIP1") ~ "Inflammasome"))
pathway_level <- c("Healthy controls","Polygenic JIA","NFKB","Inflammasome","Immune metabolism")

cell_color_vec6 <- c("Healthy controls" = "#6699CC", "Polygenic JIA" = "#FEC260",
                     "NOD2" = "#EC706E", "PSTPIP1" = "#96AF95",
                     "LPIN2" = "#FFEFC1", "LACC1" = "#C8BBD5")
cell_color_vec5 <- c("Healthy controls" = "#6699CC", "Polygenic JIA" = "#FEC260",
                     "NFKB" = "#B3D1E7", "Inflammasome" = "#B07A99",
                     "Immune metabolism" = "#847AB3")

plist_celltype <- PlotCountDistribution(pbmc1, Count_var = "celltype", color_vec = cell_color_vec)
plist_genome <- PlotCountDistribution(pbmc1, Count_var = "genome", color_vec = cell_color_vec6, group_order = genome_level)
plist_status <- PlotCountDistribution(pbmc1, Count_var = "status", color_vec = cell_color_vec3)   # original kept

plist_celltype$lollipop / plist_genome$donut
plist_celltype$donut / plist_status$donut
plist_celltype$donut / plist_genome$donut

p_nested <- PlotNestedDonut(pbmc1,
                            outer_var = "pathway", inner_var = "genome",
                            outer_color = cell_color_vec5, inner_color = cell_color_vec6,
                            outer_level = pathway_level, inner_level = genome_level)
plist_celltype$donut / p_nested

###check data
sort(table(pbmc1$celltype)/ncol(pbmc1),decreasing = T)
sort(table(pbmc1$genome)/ncol(pbmc1),decreasing = T)
sort(table(pbmc1$pathway)/ncol(pbmc1),decreasing = T)

####(2) stack proportion grouped by samples or status
library(Seurat)
library(tidyverse)
library(patchwork)

PlotStackedBar <- function(seu, Count_var, group_by_var, facet_by_status=FALSE, sample_order, color_vec){
  if(! "status" %in% colnames(seu@meta.data)) stop("metadata must contain 'status' column")
  status_level_fixed <- c("Healthy controls","Polygenic JIA","Monogenic JIA")
  if(group_by_var=="status"){
    df <- table(seu@meta.data[[Count_var]], seu@meta.data[["status"]]) %>%
      as.data.frame() %>%
      set_names(c("stack_group","x_group","n")) %>%
      group_by(x_group) %>%
      mutate(frac = n / sum(n)*100) %>%
      ungroup()
    df$x_group <- factor(df$x_group, levels = status_level_fixed)
  }else{
    df <- table(seu@meta.data[[Count_var]], seu@meta.data[[group_by_var]]) %>%
      as.data.frame() %>%
      set_names(c("stack_group","x_group","n")) %>%
      left_join(seu@meta.data[,c(group_by_var,"status")] %>% distinct(), by=c("x_group"=group_by_var)) %>%
      group_by(x_group) %>%
      mutate(frac = n / sum(n)*100) %>%
      ungroup()
    df$x_group <- factor(df$x_group, levels=sample_order)
  }
  df$stack_group <- factor(df$stack_group, levels=names(color_vec))
  if("status" %in% colnames(df)) df$status <- factor(df$status, levels=status_level_fixed)
  p_base <- ggplot(df, aes(x=x_group, y=frac, fill=stack_group)) +
    geom_col(position=position_stack(), width=0.7) +
    scale_fill_manual(values=color_vec) +
    labs(y="Percentage (%)", x=NULL) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, hjust=1, colour="black"),
          axis.text.y=element_text(colour="black"),
          legend.title=element_blank(),
          plot.margin=margin(5,5,5,5)) +
    scale_y_continuous(expand=c(0,0))
  p_no_facet <- p_base
  p_facet <- NULL
  if(isTRUE(facet_by_status) && group_by_var!="status"){
    p_facet <- p_base + facet_wrap(~status, ncol=3, scales="free_x") +
      theme(axis.text.x=element_text(angle=45, hjust=1, colour="black"))
  }
  p_horiz_no_facet <- NULL
  p_horiz_facet <- NULL
  if(group_by_var != "status"){
    df$y_group <- factor(df$x_group, levels=rev(sample_order))
    p_horiz_base <- ggplot(df, aes(x=frac, y=y_group, fill=stack_group)) +
      geom_col(position=position_stack(), width=0.7) +
      scale_fill_manual(values=color_vec) +
      labs(x="Percentage (%)", y=NULL) +
      theme_classic() +
      theme(axis.text.x=element_text(colour="black"),
            axis.text.y=element_text(colour="black"),
            legend.title=element_blank(),
            plot.margin=margin(5,5,5,5)) +
      scale_x_continuous(expand=c(0,0))
    p_horiz_no_facet <- p_horiz_base
    if(isTRUE(facet_by_status)){
      p_horiz_facet <- p_horiz_base +
        facet_wrap(~status, ncol=1, scales="free_y", strip.position="right") +
        theme(strip.placement = "outside")
    }
  }
  return(list(stack_no_facet=p_no_facet,
              stack_facet=p_facet,
              horiz_stack_no_facet=p_horiz_no_facet,
              horiz_stack_facet=p_horiz_facet))
}

# grouping columns
sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "157LACC1","182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
pbmc1$sample <- factor(pbmc1$sample, levels = sample_level)
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(genome = case_when(
  sample %in% c("C1","C2","C3","C4","C5","C6") ~ "Healthy controls",
  sample %in% c("poly1","poly2","poly3","poly4","poly5","poly6") ~ "Polygenic JIA",
  sample %in% c("157LACC1") ~ "LACC1",
  sample %in% c("182LPIN2") ~ "LPIN2",
  sample %in% c("181NOD2","183NOD2") ~ "NOD2",
  sample %in% c("190PSTPIP1","200PSTPIP1") ~ "PSTPIP1"))
genome_level <- c("Healthy controls","Polygenic JIA","NOD2","PSTPIP1","LPIN2","LACC1")
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(pathway = case_when(
  sample %in% c("C1","C2","C3","C4","C5","C6") ~ "Healthy controls",
  sample %in% c("poly1","poly2","poly3","poly4","poly5","poly6") ~ "Polygenic JIA",
  sample %in% c("157LACC1") ~ "Immune metabolism",
  sample %in% c("182LPIN2") ~ "Inflammasome",
  sample %in% c("181NOD2","183NOD2") ~ "NFKB",
  sample %in% c("190PSTPIP1","200PSTPIP1") ~ "Inflammasome"))
pathway_level <- c("Healthy controls","Polygenic JIA","NFKB","Inflammasome","Immune metabolism")
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(mutation = case_when(
  sample %in% c("C1","C2","C3","C4","C5","C6") ~ "Healthy controls",
  sample %in% c("poly1","poly2","poly3","poly4","poly5","poly6") ~ "Polygenic JIA",
  sample %in% c("181NOD2") ~ "M513K",
  sample %in% c("183NOD2") ~ "C495Y",
  sample %in% c("182LPIN2") ~ "I67R/K781TfsTer7",
  sample %in% c("190PSTPIP1","200PSTPIP1") ~ "E250K",
  sample %in% c("157LACC1") ~ "S366Ter/N19KfsTer31"))
mutation_level <- c("Healthy controls","Polygenic JIA","M513K","C495Y",
                    "I67R/K781TfsTer7","E250K","S366Ter/N19KfsTer31")
mutation_expand_level <- c("C1","C2","C3","C4","C5","C6",
                           "poly1","poly2","poly3","poly4","poly5","poly6",
                           "M513K","C495Y","I67R/K781TfsTer7","E250K","S366Ter/N19KfsTer31")
pathway_expand_level <- c("C1","C2","C3","C4","C5","C6",
                          "poly1","poly2","poly3","poly4","poly5","poly6",
                          "NFKB","Inflammasome","Immune metabolism")
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(mutation_expand = ifelse(
  sample %in% c("C1","C2","C3","C4","C5","C6","poly1","poly2","poly3","poly4","poly5","poly6"),
  as.character(sample), mutation))
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(pathway_expand = ifelse(
  sample %in% c("C1","C2","C3","C4","C5","C6","poly1","poly2","poly3","poly4","poly5","poly6"),
  as.character(sample), pathway))

# colors (cell_color_vec = your 15-color celltype palette, local)
status_level <- c("Healthy controls","Polygenic JIA","Monogenic JIA")
cell_color_vec3 <- c("#6699CC","#FEC260","#AF478A")
names(cell_color_vec3) <- status_level
cell_color_vec2 <- c("#847AB3","#96AF95","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                     "#F1BAAB","#BB0021","#D9C2D9","#FEC260","#AF478A",
                     "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9","#B1DA99",
                     "#7B1FA2","#FFEFC1")
names(cell_color_vec2) <- sample_level

##(a) Cell prop grouped by 6 HC, 6 Polygenic patients, and 6 monogenic samples
res_hf <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="sample", facet_by_status=TRUE, sample_order=sample_level, color_vec=cell_color_vec)
res_hf$horiz_stack_facet + theme(legend.position = "none")
res_hf$stack_facet + theme(legend.position = "none")

##(b) Cell prop grouped by healthy controls, polygenetic, and monogenic patients
res_status_agg <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="status", facet_by_status=FALSE, sample_order=sample_level, color_vec=cell_color_vec)
res_status_agg$stack_no_facet

##(c) Cell prop grouped by healthy controls, polygenetic, and 6 gene mutations
res_mut_big <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="mutation",
                              facet_by_status=FALSE, sample_order=mutation_level, color_vec=cell_color_vec)
res_mut_big$stack_no_facet
res_mut_big$horiz_stack_no_facet

##(d) Cell prop grouped by  6 HC, 6 Polygenic patients, and 6 gene mutations
res_mut_samp <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="mutation_expand",
                               facet_by_status=TRUE, sample_order=mutation_expand_level, color_vec=cell_color_vec)
res_mut_samp$stack_facet + theme(legend.position = "none")
res_mut_samp$horiz_stack_facet + theme(legend.position = "none")
res_mut_samp$stack_no_facet
res_mut_samp$horiz_stack_no_facet

##(e) Cell prop grouped by healthy controls, polygenetic, and  pathways
res_path_big <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="pathway",
                               facet_by_status=FALSE, sample_order=pathway_level, color_vec=cell_color_vec)
res_path_big$stack_no_facet
res_path_big$horiz_stack_no_facet

##(f) Cell prop grouped by6 HC, 6 Polygenic patients, and  pathways
res_path_samp <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="pathway_expand",
                                facet_by_status=TRUE, sample_order=pathway_expand_level, color_vec=cell_color_vec)
res_path_samp$stack_facet + theme(legend.position = "none")
res_path_samp$horiz_stack_facet + theme(legend.position = "none")
res_path_samp$stack_no_facet
res_path_samp$horiz_stack_no_facet

####(3)significant test for cell proportions
resultCount<-as.data.frame.matrix(table(pbmc1$sample,pbmc1$celltype))
resultCount$sample<-rownames(resultCount)
resultCount<-resultCount%>%select(sample,everything())
write.csv(resultCount,file="E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\resultCount.csv",row.names=F,fileEncoding = "UTF-8")
write.csv(resultCount,file = "resultCount.csv",
row.names = FALSE,fileEncoding = "UTF-8")


##(a) Compare for three groups
library(tidyverse)
sccoda_hm<-read.csv("./CellPreprocess/0827scCODA_mono_vs_HC.csv")
sccoda_hp<-read.csv("./CellPreprocess/0827scCODA_poly_vs_HC.csv")
sccoda_pm<-read.csv("./CellPreprocess/0827scCODA_mono_vs_poly.csv")
allsccoda<-do.call("rbind",list(sccoda_hm,sccoda_hp,sccoda_pm))
allsccoda<-allsccoda[allsccoda$Significant=="True",]

allsccoda$Contrast<-factor(allsccoda$Contrast,levels=c("polygenic_vs_HC","monogenic_vs_HC","monogenic_vs_polygenic"))
contrast_col<-c("polygenic_vs_HC"="#FEC260","monogenic_vs_HC"="#AF478A","monogenic_vs_polygenic"="#6699CC")

plot_scCODA_sig_lollipop_facet<-function(df_sig,point_size=2.5,errorbar_linewidth=0.35){
  df_sig$Cell.Type<-factor(df_sig$Cell.Type,levels=unique(df_sig$Cell.Type))
  
  p1<-ggplot(df_sig,aes(x=log2FC,y=Cell.Type,color=Contrast))+
    geom_vline(xintercept=0,linetype="dashed",color="black",linewidth=0.4)+
    geom_errorbarh(aes(xmin=HDI.3.,xmax=HDI.97.),height=0.22,linewidth=errorbar_linewidth)+
    geom_point(size=point_size)+
    scale_color_manual(values=contrast_col)+
    labs(x="log2FC (cell proportion)",y="Cell Type")+
    theme_classic()+
    theme(panel.border=element_rect(colour="black",fill=NA,linewidth=0.6),
          panel.grid.major=element_line(colour="grey90",linewidth=0.25,linetype="dashed"),
          panel.grid.minor=element_line(colour="grey90",linewidth=0.15,linetype="dashed"),
          axis.text.x=element_text(angle=0,hjust=0.5,colour="black"),
          axis.text.y=element_text(colour="black"),
          axis.title=element_text(colour="black"),
          axis.ticks=element_line(colour="black"),
          legend.position="none",
          strip.background=element_rect(fill="white",color="black"),
          strip.text=element_text(colour="black"))+
    facet_wrap(~Contrast,scales="free",ncol=3)
  
  p2<-ggplot(df_sig,aes(x=Cell.Type,y=log2FC,color=Contrast))+
    geom_hline(yintercept=0,linetype="dashed",color="black",linewidth=0.4)+
    geom_errorbar(aes(ymin=HDI.3.,ymax=HDI.97.),width=0.22,linewidth=errorbar_linewidth)+
    geom_point(size=point_size)+
    scale_color_manual(values=contrast_col)+
    labs(x="Cell Type",y="log2FC (cell proportion)")+
    theme_classic()+
    theme(panel.border=element_rect(colour="black",fill=NA,linewidth=0.6),
          panel.grid.major=element_line(colour="grey90",linewidth=0.25,linetype="dashed"),
          panel.grid.minor=element_line(colour="grey90",linewidth=0.15,linetype="dashed"),
          axis.text.x=element_text(angle=45,hjust=1,colour="black"),
          axis.text.y=element_text(colour="black"),
          axis.title=element_text(colour="black"),
          axis.ticks=element_line(colour="black"),
          legend.position="none",
          strip.background=element_rect(fill="white",color="black"),
          strip.text=element_text(colour="black"))+
    facet_wrap(~Contrast,scales="free",ncol=3)
  
  return(list(p_horizontal=p1,p_vertical=p2))
}

plot_list<-plot_scCODA_sig_lollipop_facet(allsccoda,point_size=1.8,errorbar_linewidth=0.5)
p_hori<-plot_list$p_horizontal
p_vert<-plot_list$p_vertical
print(p_hori)
print(p_vert)

##(b) Compare for specific monogenic groups (by pathways)
celltype_levels <- c("Naive CD4 T","Memory CD4 T",
                     "Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                     "γδT","Naive B","Memory B","Plasma","CD14 Monocyte",
                     "CD16 Monocyte","pDC","Platelet","Erythrocyte")
cell_color_vec <- c("#847AB3","#96AF95","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                    "#F1BAAB","#BB0021","#D9C2D9","#FEC260","#AF478A",
                    "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")
names(cell_color_vec) <- celltype_levels

contrast_level <- c("Polygenic JIA_vs_Healthy controls",
                    "NFKB_vs_Healthy controls",
                    "Inflammasome_vs_Healthy controls",
                    "Immune metabolism_vs_Healthy controls",
                    "NFKB_vs_Polygenic JIA",
                    "Inflammasome_vs_Polygenic JIA",
                    "Immune metabolism_vs_Polygenic JIA")
contrast_col <- c("#FEC260","#E67E22","#E74C3C","#8E44AD",
                  "#3498DB","#1ABC9C","#2E86C1")
names(contrast_col) <- contrast_level


plot_scCODA_sig_lollipop_facet <- function(df_sig, point_size = 2.5, errorbar_linewidth = 0.35){
  df_sig$Cell.Type <- factor(df_sig$Cell.Type, levels = intersect(celltype_levels, unique(df_sig$Cell.Type)))
  df_sig$Contrast <- factor(df_sig$Contrast, levels = intersect(contrast_level, unique(df_sig$Contrast)))
  p1 <- ggplot(df_sig, aes(x = log2FC, y = Cell.Type, color = Contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = HDI.3., xmax = HDI.97.), height = 0.22, linewidth = errorbar_linewidth) +
    geom_point(size = point_size) +
    scale_color_manual(values = contrast_col) +
    labs(x = "log2FC (cell proportion)", y = "Cell Type") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey90", linewidth = 0.15, linetype = "dashed"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none",
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(colour = "black")) +
    facet_wrap(~Contrast, scales = "free", ncol = 3)
  p2 <- ggplot(df_sig, aes(x = Cell.Type, y = log2FC, color = Contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = HDI.3., ymax = HDI.97.), width = 0.22, linewidth = errorbar_linewidth) +
    geom_point(size = point_size) +
    scale_color_manual(values = contrast_col) +
    labs(x = "Cell Type", y = "log2FC (cell proportion)") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey90", linewidth = 0.15, linetype = "dashed"),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none",
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(colour = "black")) +
    facet_wrap(~Contrast, scales = "free", ncol = 3)
  colby_groups_horizontal <- ggplot(df_sig, aes(x = log2FC, y = Cell.Type, color = Contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = HDI.3., xmax = HDI.97.), height = 0.22, linewidth = errorbar_linewidth,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = point_size, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = contrast_col) +
    labs(x = "log2FC (cell proportion)", y = "Cell Type") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey90", linewidth = 0.15, linetype = "dashed"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    scale_y_discrete(drop = TRUE)
  colby_groups_vertical <- ggplot(df_sig, aes(x = Cell.Type, y = log2FC, color = Contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = HDI.3., ymax = HDI.97.), width = 0.22, linewidth = errorbar_linewidth,
                  position = position_dodge(width = 0.6)) +
    geom_point(size = point_size, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = contrast_col) +
    labs(x = "Cell Type", y = "log2FC (cell proportion)") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey90", linewidth = 0.15, linetype = "dashed"),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    scale_x_discrete(drop = TRUE)
  colby_celltype_horizontal <- ggplot(df_sig, aes(x = log2FC, y = Contrast, color = Cell.Type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = HDI.3., xmax = HDI.97.), height = 0.22, linewidth = errorbar_linewidth,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = point_size, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = cell_color_vec) +
    labs(x = "log2FC (cell proportion)", y = "Contrast") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey90", linewidth = 0.15, linetype = "dashed"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    scale_y_discrete(drop = TRUE)
  colby_celltype_vertical <- ggplot(df_sig, aes(x = Contrast, y = log2FC, color = Cell.Type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = HDI.3., ymax = HDI.97.), width = 0.22, linewidth = errorbar_linewidth,
                  position = position_dodge(width = 0.6)) +
    geom_point(size = point_size, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = cell_color_vec) +
    labs(x = "Contrast", y = "log2FC (cell proportion)") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey90", linewidth = 0.15, linetype = "dashed"),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    scale_x_discrete(drop = TRUE)
  return(list(p_horizontal = p1, p_vertical = p2,
              colby_groups_horizontal = colby_groups_horizontal,
              colby_groups_vertical = colby_groups_vertical,
              colby_celltype_horizontal = colby_celltype_horizontal,
              colby_celltype_vertical = colby_celltype_vertical))
}

plot_list <- plot_scCODA_sig_lollipop_facet(allsccoda, point_size = 3, errorbar_linewidth = 0.6)
print(plot_list$p_horizontal)
print(plot_list$p_vertical)
print(plot_list$colby_groups_horizontal)
print(plot_list$colby_groups_vertical)
print(plot_list$colby_celltype_horizontal)
print(plot_list$colby_celltype_vertical)

#####7 ----Marker gene bubble-----
library(Seurat)
library(tidyverse)
library(patchwork)
PlotMultiDotPlot <- function(seu, gene_categories, color_vec){
  DefaultAssay(seu) <- "RNA"
  plot_list <- vector("list", length(gene_categories))
  for(i in seq_along(gene_categories)){
    genei <- as.character(gene_categories[[i]])
    p <- DotPlot(seu, features = genei) +
      theme(axis.line = element_blank(),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            axis.text.x = element_text(color = "black", angle = 90, face = "italic"),
            axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1, colour = "black"),
            panel.border = element_rect(color = "black", size = 1, fill = NA),
            plot.margin = margin(0,0,0,0),
            legend.position = "none",
            panel.grid.major = element_line(color = "grey80", size = 0.5, linetype = "dashed"),
            panel.grid.minor = element_line(color = "grey90", size = 0.25, linetype = "dashed")) +
      ylab(NULL) + xlab(NULL) +
      scale_color_gradientn(colors = color_vec)
    plot_list[[i]] <- p
  }
  for(i in 2:length(plot_list)){
    plot_list[[i]] <- plot_list[[i]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  out_plot <- wrap_plots(plot_list, ncol = length(plot_list)) + plot_layout(guides = "collect")
  return(list(plot_list = plot_list, combined = out_plot))
}
gene_categories <- list(
  c("CD4","SELL","CCR7","TCF7"),
  c("CD8A","CD8B","CCR7","LEF1"),
  c("CD14","S100A8","S100A9"),
  c("IGHD","TCL1A","IL4R"),
  c("NKG7","KLRF1","KLRD1"),
  c("CD8A","CD8B","GZMA","GZMB","GZMH"),
  c("AQP3","TSHZ2","CD40LG","IL7R"),
  c("TRGC1","TRGC2","CXCR6"),
  c("FOXP3","IL2RA","IKZF2"),
  c("FCRL5","FCRL2","CD24","TNFRSF13B","CD80","IGHG3"),
  c("PF4","PPBP","ITGA2B"),
  c("FCGR3A","MS4A7","LST1"),
  c("HBB","HBA2","HBA1"),
  c("PLD4","IL3RA","IRF7"),
  c("JCHAIN","XBP1","MZB1"))


library(RColorBrewer)
col2 <- colorRampPalette(c("#BEBEBE",brewer.pal(9,"YlGnBu")))(100)[1:85]
bubble_res <- PlotMultiDotPlot(pbmc1, gene_categories = gene_categories, color_vec = col2)
bubble_res$combined

#####8 ----heatmap-----
library(GSEABase)
library(pheatmap)
library(ggplotify)
library(patchwork)

####(1) Pseudobulk: sum counts, DESEQ2 normalization resembling bulk analysis
build_vst <- function(sce_obj, target_cells, de_ref_level="ctrl", zero_thresh=10, row_sum_min=10){
  
  vst_out <- list()
  
  for(ct in target_cells){
    sub <- subset(sce_obj, subset = celltype %in% ct)
    if(ncol(sub)==0){warning(ct," no cells");next}
    
    pb <- Seurat::AggregateExpression(sub, group.by="sample", assays="RNA", slot="counts")
    pb_cnt <- round(as.matrix(pb$RNA))
    
    # 2.filter genes: the number of samples with (count sum==0) >=10
    zero_per_gene <- apply(pb_cnt, 1, function(x) sum(x==0))
    pb_cnt <- pb_cnt[zero_per_gene < zero_thresh, ]
    pb_cnt <- pb_cnt[rowSums(pb_cnt)>=row_sum_min, , drop=FALSE]
    
    if(nrow(pb_cnt)<2 || ncol(pb_cnt)<2){warning(ct," matrix too small");next}
    
    # 3.构造coldata、dds，打印levels确认参考组
    samp <- colnames(pb_cnt)
    cond <- factor(ifelse(grepl("^C",samp),"ctrl","patient"), levels = c(de_ref_level,"patient"))
    coldata <- S4Vectors::DataFrame(condition=cond, row.names=samp)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(pb_cnt, coldata, design=~condition)
    cat("\n[",ct,"] condition levels: "); print(levels(dds$condition))
    
    dds <- DESeq2::DESeq(dds)
    vst_mat <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind=TRUE))
    
    vst_out[[ct]] <- vst_mat
  }
  return(vst_out)
}
vst_list <- build_vst(
  sce_obj = pbmc1,
  target_cells = unique(pbmc1$celltype),
  de_ref_level = "ctrl",
  zero_thresh = 10,
  row_sum_min = 10
)

####(2)heatmap:log-data or counts normalization
heatmap_avg_genevec <- function(sce_obj,
                                target_cells,
                                gene_vec,
                                gene_set_name,
                                plot_title_vec,
                                chvalue,
                                color_vec,
                                mode = c("avg_log","pseudobulk_vst"),
                                vst_list = NULL,
                                preprocess_mode = c("zscore_row_internal",
                                                    "zscore_row_pre",
                                                    "log2p1_zscore_row",
                                                    "log2p1_none")){
  
  mode <- match.arg(mode)
  preprocess_mode <- match.arg(preprocess_mode)
  bk <- unique(c(seq(-1.5,1.5, length=100)))
  
  if(length(target_cells)!=length(plot_title_vec)) stop("target_cells and plot_title_vec length mismatch")
  
  temp_list <- list()
  valid_ct <- c()
  
  if(mode == "avg_log"){
    for(k in seq_along(target_cells)){
      ct <- target_cells[k]
      sub_sce <- subset(sce_obj, subset = celltype %in% ct)
      if(ncol(sub_sce)==0){warning(sprintf("celltype: %s no cells, skip",ct));next}
      h_avg <- Seurat::AverageExpression(sub_sce, group.by = "sample")
      df <- as.data.frame(h_avg$RNA)
      df <- df[rownames(df) %in% gene_vec, ]
      df <- df[-which(apply(df,1,function(x){sum(x==0)})>=3), ]
      temp_list[[ct]] <- df
      valid_ct <- c(valid_ct, ct)
    }
  } else if(mode == "pseudobulk_vst"){
    if(is.null(vst_list)) stop("mode=pseudobulk_vst requires vst_list from build_vst()")
    for(ct in target_cells){
      if(!ct %in% names(vst_list)){warning(sprintf("celltype: %s not in vst_list, skip",ct));next}
      temp_list[[ct]] <- vst_list[[ct]]
      valid_ct <- c(valid_ct, ct)
    }
  }
  
  if(length(temp_list)==0) stop("no valid celltype, exit")
  gene_intersect <- Reduce(intersect, c(lapply(temp_list, rownames), list(gene_vec)))
  if(length(gene_intersect)<2) stop("too few intersect genes for heatmap")
  
  plot_list <- list()
  for(k in seq_along(target_cells)){
    ct <- target_cells[k]
    if(!ct %in% valid_ct) next
    this_title <- paste0(plot_title_vec[k]," | ",gene_set_name)
    
    if(mode == "avg_log"){
      sub_sce <- subset(sce_obj, subset = celltype %in% ct)
      h_avg <- Seurat::AverageExpression(sub_sce, group.by = "sample")
      df <- as.data.frame(h_avg$RNA)
      df <- df[gene_intersect, , drop = FALSE]
    } else {
      df <- temp_list[[ct]][gene_intersect, , drop = FALSE]
    }
    
    mat <- as.matrix(df)
    pheat_scale <- "none"
    
    if(preprocess_mode == "zscore_row_internal"){
      pheat_scale <- "row"
    }else if(preprocess_mode == "zscore_row_pre"){
      mat <- t(scale(t(mat)))
      pheat_scale <- "row"
    }else if(preprocess_mode == "log2p1_zscore_row"){
      mat <- log2(mat + 1)
      pheat_scale <- "row"
    }else if(preprocess_mode == "log2p1_none"){
      mat <- log2(mat + 1)
      pheat_scale <- "none"
    }
    
    set.seed(123)
    p <- pheatmap::pheatmap(mat,
                            scale = pheat_scale,
                            cluster_cols = FALSE,
                            cluster_rows = TRUE,
                            show_rownames = FALSE,
                            border_color = NA,
                            main = this_title,
                            fontsize = 6.5,
                            legend = TRUE,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            legend_breaks = seq(-2,2,0.5),
                            color = color_vec,
                            annotation_legend = TRUE,
                            breaks = bk,
                            cellwidth = 8,
                            cellheight = chvalue,
                            clustering_distance_rows="correlation",
                            clustering_method = "average")
    
    plot_list[[ct]] <- ggplotify::as.ggplot(p)
  }
  
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "hv")
  
  return(list(plot_list=plot_list,combined_plot=combined_plot,gene_intersect=gene_intersect,preprocess_mode_used=preprocess_mode))
}

plot_ct_all <- names(vst_list)
gene_set <- NFKB_list
gene_set_name <- "NKFB_signature"
color_custom <- colorRampPalette(c("#2166ac","white","#b2182b"))(101)
gmt_file <- "E:\\组学数据分析\\IL1R1_TLR1\\IL1R1\\HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2024.1.Hs.gmt"
gene_sets <- getGmt(gmt_file)
NFKB_list <- geneIds(gene_sets[[1]]) 
load("E://通路基因集合//gene4pathway.Rdata")
res_vst_all <- heatmap_avg_genevec(
  sce_obj = pbmc1,
  target_cells = plot_ct_all,
  gene_vec = IFN1$gene,
  gene_set_name = gene_set_name,
  plot_title_vec = plot_ct_all,
  chvalue = 2,
  color_vec = color_custom,
  mode = "pseudobulk_vst",
  vst_list = vst_list,
  preprocess_mode = "log2p1_zscore_row"
)

res_vst_all$plot_list$pDC
res_vst_all$plot_list$`CD14 Monocyte`







