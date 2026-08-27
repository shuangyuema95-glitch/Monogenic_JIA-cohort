library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(ggplot2)
library(DoubletFinder)
library(patchwork)
library(cowplot)

setwd("E:/Cohort PPT/JIA/code")
samples=c("C1","C2","C3","C4","C5","C6","182LPIN2","181NOD2","183NOD2","190PSTPIP1","200PSTPIP1")
fold2="E:/Cohort PPT/JIA/code"


samples<-c("poly1","poly2","poly3","poly4","poly5","poly6")
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
  pK <- paramSweep(s, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(pK, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  nExp <- round(0.03 * ncol(s))
  
  s <- doubletFinder(s, PCs = 1:50, pN = 0.25, pK = best.pK,
                     nExp = nExp, reuse.pANN = FALSE, sct = FALSE)
  colnames(s@meta.data)[8]<-"scDbl_class"
  s$sample<-samples[i]
  sce[[i]] <- s
  summary_df <- rbind(summary_df, c(samples[i], ncol(s)))
}

names(sce)<-samples
summary_df
saveRDS(sce,file="1_readsceList.RDS")

#####2 ----merge samples-----
sce<-readRDS("1_readsceList.RDS")
pbmc1 <- merge(x = sce[[1]], y = sce[-1], add.cell.ids = names(sce), project = "pbmc1")
table(pbmc1$sample)
table(sapply(rownames(pbmc1@meta.data), function(x){strsplit(x, split = "_")[[1]][1]}))
table(pbmc1$scDbl_class)
pbmc1 <- subset(pbmc1, scDbl_class == "Singlet")
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
length(hvgs_top1500)#1771
hvgs<-hvgs_top1500

###correct 
VariableFeatures(pbmc1) <- hvgs
length(pbmc1@assays$RNA@var.features)
g2m_genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(pbmc1))
s_genes  <- CaseMatch(search = cc.genes$s.genes,  match = rownames(pbmc1))
pbmc1 <- CellCycleScoring(pbmc1, g2m.features = g2m_genes, s.features = s_genes)
pbmc1 <- ScaleData(pbmc1, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(pbmc1), npcs = 50)
pbmc1@meta.data <- pbmc1@meta.data %>% mutate(datasets = str_split_fixed(rownames(pbmc1@meta.data), "_", 2)[,1])
pbmc1 <- RunHarmony(pbmc1, group.by.vars = "datasets", dims.use = 1:50)# theta = 2) 


#####3 ----Cell clustering-----
pbmc1 <- FindNeighbors(pbmc1, dims = 1:50, reduction = "harmony")
pbmc1 <- FindClusters(pbmc1, resolution = 0.45)
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
Count_marker()
Count_marker2()
View(sce2.markers[sce2.markers$cluster==1,])
sce2Con.marker[match(12,unique(pbmc1$seurat_clusters))]
View(sce2.markers[sce2.markers$gene%in%CD14_Monocyte,])
View(sce2.markers[sce2.markers$gene%in%CD4T,])
grep(CD4T,sce2Con.marker)


cols = c("#DCD4D5B4","red")
DefaultAssay(pbmc1)<-"RNA"
cell_markers<- list(
  Naive_CD4_T=unique(c("CD4","CCR7","SELL","LEF1","TCF7","IL7R","LTB")),
  Memory_CD4_T=unique(c("CD4","IL7R","CD44","CD27","CXCR3","RGS1","CREM","CD3D")),
  Effector_CD4_T=unique(c("CD4","GATA3","IL17A","IL17F","CXCL3","ICOS")),
  Naive_CD8_T=unique(c("CD8A","CD8B","CCR7","SELL","LEF1","TCF7","IL7R","LTB")),
  Memory_CD8_T=unique(c("CD8A","IL7R","CD27","GZMK","EOMES","CD44","KLRG1","RGS1","CREM","CD3D")),
  Effector_CD8_T=unique(c("CD8A","IFNG","GZMB","PRF1","KLRG1","NKG7","LAMP1")),
  Effective_Memory_T=unique(c("NKG7","GZMK","IFNG")),
  Treg=unique(c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18")),
  gamma_delta_T=unique(c("TRDC","TRGC1","TRGC2","RORC","ZBTB16","CXCR6")),
  NKT=unique(c("ZBTB16","TRAC","CD3D","KLRD1","KLRK1","KLRB1","IL17RA","CD4")),
  NK=unique(c("NKG7","GZMB","PRF1","KLRD1","KLRK1","KLRC1","KLRB1","CCL5")),
  B_cell=unique(c("CD19","CD79A","CD79B","MS4A1","BANK1","CD37")),
  Plasma_cell=unique(c("XBP1","MZB1","SDC1","JCHAIN","IGHA1","IGHM","IGHD")),
  CD14_Monocyte=unique(c("CD14","LYZ","LY6C2","S100A8","S100A9","VCAN","CCL2","CTSS")),
  CD16_Monocyte=unique(c("FCGR3A","MS4A7","LST1","IFITM3","CSF1R","CD68")),
  DC=unique(c("FCER1A","CD1C","XCR1")),
  Macrophage=unique(c("ADGRE1","CD68","CD63","CXCL2","CCR1","NLRP3","VEGFA")),
  Neutrophil=unique(c("FCGR3B","S100A8","S100A9","MPO","ELANE","RETNLG")),
  LDG=unique(c("FCGR3B","FUT4","CEACAM8","CSF3R","CMTM2","NCF1","SRGN","TREM1","IL1RN","IFITM2")),
  aDC=unique(c("CD40","CD80","CD86","IL12A","CCR7","LAMP3","CD83","IL1B","ITGAX")),
  pDC=unique(c("PLD4","SIGLEC6","IL3RA","LILRA4","CD123","IRF7","CLEC4C","TCF4","BST2")),
  cDC1=unique(c("XCR1","CLEC9A","BATF3","IRF8","ITGAX","CADM1")),
  cDC2=unique(c("CD209","CLEC10A","FCER1G","IRF4","KLF4","SIRPA","ITGAM")),
  mDC=unique(c("CD1C","CD11C","ITGAX","FCER1A","CLEC10A","HLA-DQB1","HLA-DRA","FCGR2B","VCAN")),
  LAMP3_matureDC=unique(c("LAMP3","CCR7","FSCN1","CCL19","CCL22")),
  Platelet=unique(c("PF4","PPBP","ITGA2B","GP1BA","TUBB1")),
  CD8T=unique(c("CD3D","CD3G","CD3E","CD247","CD8A","CD8B","GZMK")),
  CD4T=unique(c("CD3D","CD3G","CD3E","CD247","CCR7","IL7R","CD4")),
  Erythrocyte=unique(c("HBB","HBA1","HBA2","ALAS2","ANK1")),
  Cytotoxic_CD8T=unique(c( "CD8A","CD8B","NKG7","GNLY","GZMB","GZMH","GZMK",
                           "PRF1","CTSW","IFNG","CCL5","TBX21","EOMES","CX3CR1","FGFBP2","LAMP1","CST7")))

B_markers <- list(
  # Pan-B cell
  Total_B = c("CD19", "CD79A", "CD79B", "MS4A1", "CD37", "CD74", "HLA-DRA"),
  # Transitional B: CD24hi CD38hi CD10+ IgMhi CD27-
  Transitional_B = c("TCL1A", "CD24", "CD38", "MME", "IL4R", "IGHM", "IGHD"),
  # Mature naive B: IgM+ IgDhi CD27- CD38low
  Naive_B = c("IGHM", "IGHD", "FCER2", "MS4A1", "CD22", "HVCN1", "TCL1A"),
  # Unswitched memory B: CD27+ IgM+ ±IgD
  Unswitched_memory_B = c("CD27", "IGHM", "IGHD", "FCRL1", "TNFRSF13B", "GPR183"),
  # Switched memory B: IgD- IgG+ or IgA+, usually CD27+
  Switched_memory_B = c("CD27", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "TNFRSF13B", "GPR183"),
  # Double-negative B: CD27- IgD-
  Double_negative_B = c("FCRL5", "FCRL3", "ITGAX", "TBX21", "IGHG1", "IGHA1"),
  # Atypical memory / DN2-like B
  Atypical_memory_B = c("FCRL5", "FCRL3", "ITGAX", "TBX21", "SIGLEC6", "ZEB2"),
  # Activated B: activation/state markers
  Activated_B = c("CD69", "CD83", "CD86", "CD40", "MYC", "NFKBIA", "TNFAIP3"),
  # Plasmablast: antibody-secreting program
  Plasmablast = c("CD38", "MZB1", "JCHAIN", "XBP1", "PRDM1", "IRF4", "DERL3", "FKBP11", "SEC11C"),
  # Mature plasma cell: terminally differentiated ASC
  Plasma_cell = c("SDC1", "TNFRSF17", "MZB1", "JCHAIN", "XBP1", "PRDM1", "IRF4", "DERL3", "FKBP11", "SEC11C"))
DC_markers <- list(
  mDC_cDC = c("ITGAX","CD1C","HLA-DRA","CD86","CD83"),
  cDC1 = c("ITGAX","XCR1","CLEC9A"),
  cDC2 = c("ITGAX","CD1C","SIRPA"),
  pDC = c("TCF4","IRF7","CLEC4C","LILRA4"))

DefaultAssay(pbmc1)<-"RNA"
DotPlot(pbmc1,features=unique(
  c(DC_markers$cDC2
    )),
  assay='RNA',
  cols=cols)+RotatedAxis()+theme(
    axis.line = element_blank(),
    legend.text = element_text(size=10),
    legend.title=element_text(size=10),
    axis.text.x = element_text(angle= 90 , vjust= .5 , hjust= 1 ),
    axis.text.y= element_text(angle= 45 , vjust= .5 , hjust= 1 ),
    panel.border = element_rect(color = "black", size = 1, fill = NA))+ylab(NULL)+xlab(NULL)


Idents(pbmc1) <- "seurat_clusters"
pbmc1<-RenameIdents(pbmc1,
                    '0'='Naive CD4 T',#"CD4","SELL","CCR7","LEF1" done 2
                    '1'='Naive CD8 T',#"CD8A","CD8B","CCR7","LEF1" done 2
                    '2'='CD14 Monocyte',#"CD14","S100A8","S100A9" done 2
                    '3'='Naive B',#"CD19","IGHM","IGHD","TCL1A" done 2
                    '4'='NK', #"NKG7","gzmb","KLRD1" done   2       
                    '5'='Cytotoxic CD8 T',# "CD8A","CD8B","GZMK","GZMH" done 2
                    '6'='Naive CD4 T',#"CD4","SELL","CCR7","LEF1" done 2
                    '7' = 'Naive CD8 T',#"CD8A","CD8B","CCR7","LEF1" done 2
                    '8' = 'Treg',#"FOXP3","IL2RA","IKZF2" done 2
                    '9'='Platelet',#"PF4","PPBP","ITGA2B" done 2
                    '10'='Memory B',# "GPR183","TNFRSF13B","CD24" done 2
                    '11'='CD14 Monocyte',#2
                    '12'='γδT',#"TRGC1", "TRGC2","CXCR6" done 2
                    '13'='Erythrocyte', # "HBB","HBA1","HBA2" done 2         
                    '14'='CD16 Monocyte',#"FCGR3A","MS4A7","LST1" done 2
                     '15'='Naive B',#"CD19","IGHM","IGHD","TCL1A" done 2
                    '16'='NK',#"NKG7","gzmb","KLRD1" done   2 
                    '17'='CD14 Monocyte',#"CD14","S100A8","S100A9" done 2
                    '18'='pDC',#"PLD4","IL3RA","IRF7" done 2
                    '19'='Plasma',#"JCHAIN","XBP1","MZB1" CD19- CD38+ CD27+ done 2
                    '20'='CD14 Monocyte',#"CD14","S100A8","S100A9" done 2
                    '21'='Naive B'#"CD19","IGHM","IGHD","TCL1A" done 2
                    )
pbmc1@meta.data$celltype <- Idents(pbmc1)
length(unique(pbmc1$celltype))#15
sample_map <- tibble(
  sample = c("C1","C2","C3","C4","C5","C6",
             "poly1","poly2","poly3","poly4","poly5","poly6",
             "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1"),
  sample2 = c("C1","C2","C3","C4","C5","C6",
              "PJ1","PJ2","PJ3","PJ4","PJ5","PJ6",
              "P9","P12","P13","P14","P15"))
pbmc1@meta.data <- pbmc1@meta.data %>%left_join(sample_map, by = "sample")

pbmc1@meta.data<-pbmc1@meta.data%>%mutate(status=case_when(
  sample%in%c("C1","C2","C3","C4","C5","C6")~"Healthy controls",
  sample%in%c("poly1","poly2","poly3","poly4","poly5","poly6")~"Polygenic JIA",
  sample%in%c("182LPIN2","181NOD2","183NOD2","190PSTPIP1","200PSTPIP1")~"Monogenic JIA"))
status_level<-c("Healthy controls","Polygenic JIA","Monogenic JIA")

#####5 ----Cell UMAP-----
celltype_levels <- c("Naive CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                     "γδT","Naive B","Memory B","Plasma","CD14 Monocyte",
                     "CD16 Monocyte","pDC","Platelet","Erythrocyte")
cell_color_vec <- c("#847AB3","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                    "#F1BAAB","#BB0021", "#D9C2D9","#FEC260","#AF478A",
                    "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")
names(cell_color_vec) <- celltype_levels

sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
cell_color_vec2 <- c("#847AB3","#6699CC","#86C7B4","#007ABA","#9CD2ED",
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
ggsave(p1,file="E:/Cohort PPT/JIA/code//00_0.05umap.pdf",width = 19.97,height = 3.25)



#####6 ----Cell Proportion-----
####(1) total cell count to calculated cell fraction
PlotCountDistribution <- function(seu, Count_var, color_vec){
  df <- table(seu@meta.data[[Count_var]]) %>%
    as.data.frame() %>%
    set_names(c("group","n")) %>%
    mutate(frac = n / sum(n) * 100, label = sprintf("%.1f%%", frac)) %>%
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
  left = 0.32, bottom = 0.018, right = 0.98, top = 0.38
)
p_comb 
ggsave(p_comb,file="E:/Cohort PPT/JIA/code//01_PropDonut.pdf",width = 3.27,height = 3.38)

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
  
  # original vertical stacked
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
  
  # horizontal flip: x=Percentage, y=sample; sample top‑down = sample_order
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

sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")
cell_color_vec2 <- c("#847AB3","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                     "#F1BAAB","#BB0021", "#D9C2D9","#FEC260","#AF478A",
                     "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9","#B1DA99",
                     "#7B1FA2","#FFEFC1")
names(cell_color_vec2) <- sample_level
pbmc1$sample <- factor(pbmc1$sample, levels = sample_level)

status_level <- c("Healthy controls","Polygenic JIA","Monogenic JIA")
cell_color_vec3 <- c("#6699CC","#FEC260","#AF478A")
names(cell_color_vec3) <- status_level

res_hf <- PlotStackedBar(pbmc1, Count_var="celltype", group_by_var="sample", facet_by_status=TRUE, sample_order=sample_level, color_vec=cell_color_vec)
res_hf$horiz_stack_facet+theme(legend.position = "none")
res_hf$stack_facet+theme(legend.position = "none")

res_status_agg <- PlotStackedBar(pbmc1,Count_var = "celltype",group_by_var = "status",facet_by_status = FALSE,sample_order = sample_level,color_vec = cell_color_vec)
res_status_agg$stack_no_facet

####(3)significant test for cell proportions
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
  c("CD4","SELL","CCR7","LEF1"),
  c("CD8A","CD8B"),
  c("CD14","S100A8","S100A9"),
  c("TCL1A","CD19","IGHM","IGHD"),
  c("NKG7","GZMB","KLRD1"),
  c("CD8A","CD8B","GZMK","GZMH"),
  c("FOXP3","IL2RA","IKZF2"),
  c("PF4","PPBP","ITGA2B"),
  c("GPR183","TNFRSF13B","CD24"), 
  c("TRGC1", "TRGC2","CXCR6"),
  c("HBB","HBA1","HBA2"),
  c("FCGR3A","MS4A7","LST1"),
  c("PLD4","IL3RA","IRF7"),
  c("JCHAIN","XBP1","MZB1"))
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







