library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ggplotify)
library(cowplot)
library(patchwork)

sample_level <- c("C1","C2","C3","C4","C5","C6",
                  "poly1","poly2","poly3","poly4","poly5","poly6",
                  "182LPIN2","181NOD2","183NOD2","200PSTPIP1","190PSTPIP1")

mono_metab_sub <- subset(pbmc1, celltype %in% c("CD14 Monocyte","CD16 Monocyte"))
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

keep_pathways<-c("Glycolysis / Gluconeogenesis","Oxidative phosphorylation","Fatty acid biosynthesis","Fatty acid degradation","Glycerolipid metabolism","Glycerophospholipid metabolism","Ether lipid metabolism","Sphingolipid metabolism","Arachidonic acid metabolism","Pyruvate metabolism","Citrate cycle (TCA cycle)","Pentose phosphate pathway","Glutathione metabolism")

plot_metab_heatmap<-function(data,cells,col_palette,bk_vals,cellw,cellh,keep_pathways,sample_level){
  allheat<-list()
  for(i in seq_along(cells)){
    subdata<-data%>%filter(celltype==cells[i])
    avg_data<-subdata%>%pivot_longer(cols=all_of(keep_pathways),names_to="Pathway",values_to="Score")%>%
      group_by(sample,Pathway)%>%summarise(mean_score=mean(Score,na.rm=TRUE),.groups="drop")%>%
      pivot_wider(names_from=sample,values_from=mean_score)
    mat<-avg_data%>%column_to_rownames("Pathway")%>%as.matrix()
    mat<-mat[,intersect(sample_level,colnames(mat)),drop=FALSE]
    allheat[[i]]<-as.ggplot(pheatmap(mat,scale="row",cluster_cols=F,cluster_rows=T,show_rownames=T,
                                     border_color=NA,main=cells[i],fontsize=6,legend=T,treeheight_col=0,treeheight_row=0,
                                     legend_breaks=seq(-2,2,0.5),color=col_palette,annotation_legend=T,breaks=bk_vals,
                                     cellwidth=cellw,cellheight=cellh,fontsize_col=8,fontsize_row=10,
                                     clustering_distance_rows="euclidean",clustering_method="median"))
  }
  return(allheat)
}

# pre‑filter input df: retain celltype,sample and keep_pathways only
mat<-metab_score_df[,intersect(c("celltype","sample",keep_pathways),colnames(metab_score_df))]
bk_metab<-unique(c(seq(-1.5,1.5,length=100)))
cells<-c("CD14 Monocyte","CD16 Monocyte")

col2<-colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
# col2<-heat_col <- colorRampPalette(c("#007ABA", "#FFFFFF", "#AF478A"))(101)#"#007ABA""#6699CC"
# col2<-heat_col <- colorRampPalette(c("#007ABA", "#FFFFFF", "#BB0021"))(101)
# col2<- colorRampPalette(c("#847AB3", "#FFFFFF", "#EC706E"))(101)

Metab_heats<-plot_metab_heatmap(mat,cells,col_palette=col2,
                                bk_vals=bk_metab,
                                cellw=16,cellh=8.35,
                                keep_pathways=keep_pathways,sample_level=sample_level)
CD14_mono_heat<-Metab_heats[[1]]
CD16_mono_heat<-Metab_heats[[2]]
CD14_mono_heat/CD16_mono_heat
