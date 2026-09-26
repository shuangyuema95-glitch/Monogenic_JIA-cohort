library(ggplot2)
library(openxlsx)
library(tidyverse)


clinidata<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
grep("ystem",unique(clinidata$Number),value = T)
unique(as.character(as.matrix(clinidata[13:92,3:ncol(clinidata)])))#92-13+1=80

######1 system prop
sysdata<-clinidata[clinidata$Number%in%grep("ystem",unique(clinidata$Number),value = T),]
dim(sysdata)#80 35
colnames(sysdata)
exist<-data.frame()
alle<-data.frame()
pt_cols <- grep("^P\\d+$", colnames(sysdata), value = TRUE)
for(pt in pt_cols){
  allP<-sysdata[,c("Number",pt)]
  yes_sys<-unique(allP[allP[,2]=="yes",'Number'])
  eval_sys<-unique(allP[allP[,2]!="not available",'Number'])
  sys<-data.frame(value=unique(grep("ystem",unique(clinidata$Number),value = T)))
  sys$count<-ifelse(sys$value%in%yes_sys,1,0)
  sys1<-as.data.frame((t(sys)))
  colnames(sys1)<-sys1['value',]
  sys1<-sys1[-1,]
  exist<-rbind(exist,sys1)
  sys2<-data.frame(value=unique(grep("ystem",unique(clinidata$Number),value = T)))
  sys2$count<-ifelse(sys2$value%in%eval_sys,1,0)
  sys2<-as.data.frame((t(sys2)))
  colnames(sys2)<-sys2['value',]
  sys2<-sys2[-1,]
  alle<-rbind(alle,sys2)
}
exist1<-as.data.frame(apply(exist,2,as.numeric))
alle1<-as.data.frame(apply(alle,2,as.numeric))
rownames(exist1)<-pt_cols
rownames(alle1)<-pt_cols
pheno_prop <- data.frame(
  System = colnames(exist1),
  Positive_n = colSums(exist1),
  Evaluated_n = colSums(alle1),
  Percentage = round(colSums(exist1)/colSums(alle1)*100, 2)
)
pheno_prop

#####(1) check again for this analysis, Clinical system involvement (IBD list style)
library(readxl)
clinidata <- read.xlsx("E:/Cohort PPT/JIA/code/clinical data.xlsx")
pt_cols <- grep("^P\\d+$", colnames(clinidata), value = TRUE)
sys_list <- grep("ystem", unique(clinidata$Number), value = TRUE)
clin <- clinidata[clinidata$Number %in% sys_list, c("Number","Number_ID",pt_cols)]
clin$System <- NA
cur <- NA
for(k in 1:nrow(clin)){
  if(clin$Number[k] %in% sys_list) cur <- clin$Number[k]
  clin$System[k] <- cur
}
system_map <- list()
for(s in sys_list){
  system_map[[s]] <- clin$Number_ID[clin$System == s]
}
mat <- as.matrix(clin[, pt_cols])
rownames(mat) <- clin$Number_ID
T_pM <- as.data.frame(t(mat), stringsAsFactors = FALSE)
reSys <- list(); allSys <- list()
for(i in 1:nrow(T_pM)){
  allP <- unlist(T_pM[i, ])
  yesP <- names(allP[allP == "yes"])
  res <- unlist(lapply(system_map, function(x) length(intersect(x, yesP))))
  res[res != 0] <- 1
  reSys[[i]] <- res
  al <- sapply(system_map, function(phe_set){
    any(allP[phe_set] != "not available")
  })
  allSys[[i]] <- as.integer(al)
}
reSys <- do.call("rbind", reSys)
allSys <- do.call("rbind", allSys)
colnames(allSys) <- names(system_map)
pheno_prop <- data.frame(
  System = names(colSums(reSys)),
  Positive_n = colSums(reSys),
  Evaluated_n = colSums(allSys),
  Percentage = round(colSums(reSys)/colSums(allSys)*100, 2)
)
pheno_prop

#####(2) total prop bar
pheno_prop <- pheno_prop[order(-pheno_prop$Percentage), ]
pheno_prop$System <- factor(pheno_prop$System, levels = rev(pheno_prop$System))
prop_bar<-ggplot(pheno_prop, aes(x = System, y = Percentage)) +
  geom_col(fill = "#847AB3", width = 0.7) +
  geom_text(aes(label = paste0(Percentage, "%"), y = Percentage + 2), size = 3, color = "black") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  labs(x = "", y = "Percentage (%)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 8))

prop_bar

#####(3)Pie chart for total proportion
pheno_prop <- pheno_prop[order(-pheno_prop$Percentage), ]
color_vec <- c("#847AB3","#6699CC","#86C7B4","#007ABA","#9CD2ED",
               "#F1BAAB","#BB0021","#D9C2D9","#FEC260","#AF478A",
               "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")
pheno_prop$color <- color_vec[1:nrow(pheno_prop)]
pheno_prop$label <- paste0(pheno_prop$System, " (", pheno_prop$Percentage, "%)")
pheno_prop$System <- factor(pheno_prop$System, levels = pheno_prop$System)
ggplot(pheno_prop, aes(x = "", y = Percentage, fill = System)) +
  geom_col(width = 1, color = NA) +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(Percentage, "%"),
                y = Percentage/2 + cumsum(c(0, head(Percentage, -1)))),
            size = 3, color = "white") +
  scale_fill_manual(values = setNames(pheno_prop$color, pheno_prop$System),
                    breaks = pheno_prop$System,
                    labels = pheno_prop$label) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8, color = "black"))



######2 system prop group by  pathway
clinidata<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
metainfo<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")
as.character(as.matrix(clinidata[2,3:35]))->x1
metainfo$gene->x2
identical(x1,x2)

as.data.frame(t(data.frame(c("category","pathway",metainfo$pathway))))->y1
colnames(y1)<-colnames(clinidata)
clinidata2<-rbind(clinidata,y1)
rownames(clinidata2)[93]<-93

Count_prop<-function(sysdata){
  exist<-data.frame()
  alle<-data.frame()
  pt_cols <- grep("^P\\d+$", colnames(sysdata), value = TRUE)
  for(pt in pt_cols){
    allP<-sysdata[,c("Number",pt)]
    yes_sys<-unique(allP[allP[,2]=="yes",'Number'])
    eval_sys<-unique(allP[allP[,2]!="not available",'Number'])
    sys<-data.frame(value=unique(grep("ystem",unique(clinidata$Number),value = T)))
    sys$count<-ifelse(sys$value%in%yes_sys,1,0)
    sys1<-as.data.frame((t(sys)))
    colnames(sys1)<-sys1['value',]
    sys1<-sys1[-1,]
    exist<-rbind(exist,sys1)
    sys2<-data.frame(value=unique(grep("ystem",unique(clinidata$Number),value = T)))
    sys2$count<-ifelse(sys2$value%in%eval_sys,1,0)
    sys2<-as.data.frame((t(sys2)))
    colnames(sys2)<-sys2['value',]
    sys2<-sys2[-1,]
    alle<-rbind(alle,sys2)
  }
  exist1<-as.data.frame(apply(exist,2,as.numeric))
  alle1<-as.data.frame(apply(alle,2,as.numeric))
  rownames(exist1)<-pt_cols
  rownames(alle1)<-pt_cols
  pheno_prop <- data.frame(
    System = colnames(exist1),
    Positive_n = colSums(exist1),
    Evaluated_n = colSums(alle1),
    Percentage = round(colSums(exist1)/colSums(alle1)*100, 2)
  )
  return(pheno_prop)
}
pathways<-unique(as.character(as.matrix(clinidata2[93,3:35])))

pathway_sysprop<-data.frame()
for(i in 1:length(pathways)){
  index1<-pathways[i]
  grep(index1,clinidata2[93,])->colindex
  sysdata1<-clinidata2[clinidata2$Number%in%grep("ystem",unique(clinidata2$Number),value = T),
                       c(1,2,colindex)]
  path1prop<-Count_prop(sysdata1)
  path1prop$pathway<-index1
  
  pathway_sysprop<-rbind(pathway_sysprop,path1prop)
}
View(pathway_sysprop)
#####(1) check again for this analysis,  by pathway
pathway_row <- clinidata2[93, ]
pt_cols <- grep("^P\\d+$", colnames(clinidata2), value = TRUE)
pt_pathway <- setNames(unlist(pathway_row[, pt_cols]), pt_cols)
pathways <- unique(pt_pathway)
sys_list <- grep("ystem", unique(clinidata2$Number), value = TRUE)
pheno_by_pathway <- data.frame()
for(pw in pathways){
  pts_pw <- names(pt_pathway[pt_pathway == pw])
  exist_pw <- exist1[pts_pw, , drop = FALSE]
  alle_pw <- alle1[pts_pw, , drop = FALSE]
  tmp <- data.frame(
    Pathway = pw,
    System = colnames(exist_pw),
    Positive_n = colSums(exist_pw),
    Evaluated_n = colSums(alle_pw),
    n_patients = length(pts_pw)
  )
  tmp$Percentage <- round(tmp$Positive_n / tmp$Evaluated_n * 100, 2)
  pheno_by_pathway <- rbind(pheno_by_pathway, tmp)
}
pheno_by_pathway

write.table(pheno_by_pathway, "E:/Cohort PPT/JIA/code/CellPreprocess/clinical_3d.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

#####(2) go to python and plot 3D bar
top3 <- pheno_by_pathway %>%
  group_by(Pathway) %>%
  arrange(desc(Percentage), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  as.data.frame()
top3

##save result
save(pheno_by_pathway,pheno_prop,file="E://Cohort PPT//JIA//code//Phenotype_prop.Rdata")


######3 feature proportion for each feature
sysorder<-c("Growth and musculoskeletal system","Systemic Manifestations","Hematologic system","Immune system","Cutaneous-Mucosal system","Digestive system","Hepatosplenic system","Endocrine system","Respiratory system","Nervous system","Cardiovascular system","Urinary system")
clinidata<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
grep("ystem",unique(clinidata$Number),value = T)
unique(as.character(as.matrix(clinidata[13:92,3:ncol(clinidata)])))#92-13+1=80
sysdata<-clinidata[clinidata$Number%in%grep("ystem",unique(clinidata$Number),value = T),]
dim(sysdata)#79 35
colnames(sysdata)

features_prop<-data.frame()
for(i in 1:nrow(sysdata)){
  features<-as.character(as.matrix(sysdata[i,3:ncol(sysdata)]))
  total<-sum(features!="not available")
  present<-sum(features=="yes")
  ratio<-present/total
  features_prop<-rbind(features_prop,c(as.character(as.matrix(sysdata[i,1:2])),round(ratio,2)))
}
colnames(features_prop)<-c("sys","fea","ratio")
features_prop$ratio<-round(as.numeric(features_prop$ratio),2)

as.numeric(which(apply(sysdata,1,function(x){sum(x[3:length(x)]=="yes")})==0))->index
sysdata$Number_ID[index]#26 
dim(features_prop[features_prop$ratio==0,])#26 

features_prop$fea <- trimws(features_prop$fea)
features_prop$fea <- gsub("’", "'", features_prop$fea)
features_prop$fea <- sapply(features_prop$fea, function(x) {
  s <- tolower(x)
  substr(s, 1, 1) <- toupper(substr(s, 1, 1))
  s
})
features_prop$sys <- factor(features_prop$sys, levels = sysorder)
features_prop <- features_prop[order(features_prop$sys, -features_prop$ratio), ]
features_prop$fea <- factor(features_prop$fea, levels = unique(features_prop$fea[order(features_prop$sys, -features_prop$ratio)]))

feature_prop <- ggplot(features_prop[features_prop$ratio!=0,], aes(x = fea, y = ratio)) +
  geom_col(fill = "#847AB3", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  facet_grid(. ~ sys, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 8)) +
  labs(x = "", y = "Proportion (%)")
feature_prop
grep("/", features_prop$fea, value = TRUE)

######4 gene-phenotype heatmap,grouped by genes
clinidata<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
metainfo<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")
grep("ystem",unique(clinidata$Number),value = T)
unique(as.character(as.matrix(clinidata[13:92,3:ncol(clinidata)])))#92-13+1=80
sysdata<-clinidata[clinidata$Number%in%grep("ystem",unique(clinidata$Number),value = T),]
dim(sysdata)#80 35

gene <- c("ADA2","LACC1","ELF4","XIAP","MVK","LPIN2","NLRP3","NOD2",
                "PSTPIP1","RAC2","COPA","UNC93B1","SOCS1","TNFAIP3")
pathways<-c("NFKB","Inflammasome","Uncategoried","Immune metabolism",
            "Cell death","Interferon")

L0<-list()
for(i in 1:length(pathways)){
  p<-pathways[i];filter(metainfo,pathway==p)%>%pull(gene)%>%unique()->g
 
  L1<-list()
  for(j in 1:length(g)){
    index_g<-g[j]
    pats<-metainfo[metainfo$gene%in%index_g,'ID']
    data1<-sysdata[,c(1,2,match(pats,colnames(sysdata)))]
    data1$result<-apply(data1,1,function(x){
      y<-as.character(x)[-c(1,2)]
      any(y=="yes")->t1
      any(y=="no")->t2
      all(y=="not available")->t3
      if(t1){
        res<-c(1)
      }else if(t2){
        res<-c(0)
      }else{
        res<-c(NA)
      }
      return(res)
    })
    data1[,c("Number","Number_ID","result")]->data2
    data2$gene<-index_g
    L1[[j]]<-data2
  }
  l2<-do.call("rbind",L1)
  L0[[i]]<-l2
}
syslist<-do.call("rbind",L0)


plot_pheno_heatmap <- function(syslist, metainfo, sysdata,
                               sys_order, pw_order, ann_colors,
                               color_present = "#BB0021", color_absent = "#F7F7F7",
                               color_na = "grey80", border_color = "black",
                               cellheight = 12, cellwidth = 12){
  library(pheatmap); library(dplyr); library(tidyr); library(ggplot2)
  library(grid); library(gtable)
  syslist$result <- as.numeric(syslist$result)
  mat <- syslist %>% select(gene, Number_ID, result) %>%
    pivot_wider(names_from = Number_ID, values_from = result) %>%
    tibble::column_to_rownames("gene") %>% as.matrix()
  row_anno <- metainfo %>% select(gene, inheritance, exp, pathway, mechanism) %>%
    distinct(gene, .keep_all = TRUE) %>% tibble::column_to_rownames("gene")
  col_anno <- sysdata[, c("Number_ID", "Number")] %>% distinct(Number_ID, .keep_all = TRUE) %>%
    tibble::column_to_rownames("Number_ID")
  row_anno$pathway <- factor(row_anno$pathway, levels = pw_order)
  row_anno <- row_anno[order(row_anno$pathway), ]
  mat <- mat[rownames(row_anno), ]
  keep <- colSums(mat == 1, na.rm = TRUE) > 0
  mat <- mat[, keep]; col_anno <- col_anno[colnames(mat), , drop = FALSE]
  pw_counts <- table(row_anno$pathway)
  gaps_row <- cumsum(pw_counts[pw_counts > 0]); gaps_row <- gaps_row[-length(gaps_row)]
  col_anno$Number <- factor(col_anno$Number, levels = sys_order)
  col_anno <- col_anno[order(col_anno$Number), , drop = FALSE]
  mat <- mat[, rownames(col_anno)]
  gaps_col <- cumsum(table(col_anno$Number)[table(col_anno$Number) > 0])
  gaps_col <- gaps_col[-length(gaps_col)]
  remove_anno_border <- function(p){
    lay <- p$gtable$layout
    for(i in seq_len(nrow(lay))){
      if(grepl("annotation", lay$name[i])){
        g <- p$gtable$grobs[[i]]
        if(inherits(g, "rect")){ g$gp$col <- NA; p$gtable$grobs[[i]] <- g }
        if(!is.null(g$children)){
          for(k in seq_along(g$children)){
            if(inherits(g$children[[k]], "rect")) g$children[[k]]$gp$col <- NA
          }
        }
      }
    }
    return(p)
  }
  shrink_legend <- function(p){
    idx <- grep("annotation_legend", p$gtable$layout$name)
    for(i in idx){
      g <- p$gtable$grobs[[i]]
      if(inherits(g, "gtable")){
        g$widths[1] <- unit(0.4, "cm")
        if(length(g$widths) > 2) g$widths[2:(length(g$widths)-1)] <- unit(0.05, "cm")
        for(k in seq_along(g$grobs)){
          if(inherits(g$grobs[[k]], "rect")){
            g$grobs[[k]]$width <- unit(0.4, "cm")
            g$grobs[[k]]$height <- unit(0.4, "cm")
          }
        }
      }
      p$gtable$grobs[[i]] <- g
    }
    return(p)
  }
  common <- list(color = colorRampPalette(c(color_absent, color_present))(100),
                 annotation_row = row_anno, annotation_col = col_anno,
                 annotation_colors = ann_colors, border_color = border_color,
                 na_col = color_na, fontsize_row = 9, fontsize_col = 7,
                 cellheight = cellheight, cellwidth = cellwidth, legend = FALSE)
  p1 <- do.call(pheatmap, c(list(mat = mat, cluster_cols = TRUE, cluster_rows = FALSE,
                                 gaps_row = gaps_row, main = NA), common))
  p2 <- do.call(pheatmap, c(list(mat = mat, cluster_cols = FALSE, cluster_rows = FALSE,
                                 gaps_row = gaps_row, gaps_col = gaps_col, main = NA), common))
  p3 <- do.call(pheatmap, c(list(mat = mat, cluster_cols = FALSE, cluster_rows = TRUE,
                                 gaps_col = gaps_col, main = NA), common))
  p1 <- remove_anno_border(p1); p2 <- remove_anno_border(p2); p3 <- remove_anno_border(p3)
  p1 <- shrink_legend(p1); p2 <- shrink_legend(p2); p3 <- shrink_legend(p3)
  sys_long <- as.data.frame(mat) %>% tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "Number_ID", values_to = "result") %>%
    filter(result == 1) %>%
    left_join(rownames_to_column(col_anno, "Number_ID"), by = "Number_ID")
  stack_dat <- sys_long %>% count(gene, Number) %>%
    group_by(gene) %>% mutate(total = sum(n), prop = n/total) %>% ungroup()
  stack_dat <- stack_dat %>% group_by(gene) %>%
    mutate(rk = rank(-prop, ties.method = "first")) %>% ungroup()
  stack_dat$gene <- factor(stack_dat$gene, levels = rownames(mat))
  stack_dat$Number <- factor(stack_dat$Number, levels = sys_order)
  p_stack <- ggplot(stack_dat, aes(x = prop, y = gene, fill = Number)) +
    geom_col(width = 0.9, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(rk <= 3, round(prop*100, 2), "")),
              position = position_stack(vjust = 0.5), color = "black", size = 1.5) +
    scale_fill_manual(values = ann_colors$Number) +
    scale_x_continuous(expand = c(0,0), labels = scales::percent) +
    scale_y_discrete(expand = c(0,0)) +
    labs(x = "Proportion") +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          legend.position = "none", plot.margin = margin(0,0,0,0))
  gt <- p2$gtable
  mat_pos <- gt$layout[grep("matrix", gt$layout$name)[1], ]
  gt <- gtable_add_cols(gt, unit(3, "cm"), pos = mat_pos$r)
  gt <- gtable_add_grob(gt, ggplotGrob(p_stack), t = mat_pos$t, b = mat_pos$b,
                        l = mat_pos$r + 1, r = mat_pos$r + 1, name = "stack")
  p2$gtable <- gt
  return(list(cluster_col = p1, block = p2, cluster_row = p3, prop_table = stack_dat))
}
sys_order <- c("Growth and musculoskeletal system","Systemic Manifestations","Immune system",
               "Cutaneous-Mucosal system","Digestive system","Hepatosplenic system",
               "Endocrine system","Hematologic system","Respiratory system",
               "Nervous system","Cardiovascular system","Urinary system")
pw_order <- c("NFKB","Inflammasome","Uncategoried","Immune metabolism","Cell death","Interferon")
ann_colors <- list(
  pathway = c("Interferon"="#EC706E","NFKB"="#B3D1E7","Inflammasome"="#6699CC",
              "Immune metabolism"="#847AB3","Cell death"="#FEC260","Uncategoried"="#EBB1A4"),
  inheritance = c("AR"="#6699CC","AD"="#FEC260","XLR"="#AF478A"),
  exp = c("Balanced"="#7E75AB","Myeloid"="#82C0AE"),
  mechanism = c("deficiency"="#007ABA","gain-of-function"="#EC706E",
                "haploinsufficiency"="#F1BAAB","NA"="grey80"),
  Number = c("Growth and musculoskeletal system"="#FEC260",
             "Systemic Manifestations"="#0B71AB",
             "Immune system"="#847AB3",
             "Cutaneous-Mucosal system"="#EC706E",
             "Digestive system"="#86C7B4",
             "Hepatosplenic system"="#F1BAAB",
             "Endocrine system"="#D9C2D9",
             "Hematologic system"="#9CD2ED",
             "Respiratory system"="#FFEFC1",
             "Nervous system"="#7F7FD5",
             "Cardiovascular system"="#AF478A",
             "Urinary system"="#CBDAA9"))
res <- plot_pheno_heatmap(syslist, metainfo, sysdata, sys_order, pw_order, ann_colors,
                          color_present = "#BB0021", color_absent = "gray90",
                          color_na = "grey80", border_color = "white",
                          cellheight = 12, cellwidth = 12)
#res$cluster_col
res$block
#res$cluster_row
#res$prop_table
library(ggplotify)
heat1<-as.ggplot(res$block)
heat1
