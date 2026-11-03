library(ggplot2)
library(openxlsx)
library(tidyverse)
metainfo<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx");dim(metainfo)
country_df <- metainfo %>%
  count(country, name = "n") %>%
  mutate(percent = n / sum(n) * 100)

#####1 Figure 1  country proportion pie
country_colors <- c(
  "China"  = "gray90",
  "Egypt"  = "#B5BBE3",
  "Brazil" = "#CBDAA9")
p_country <- ggplot(country_df, aes(x = "", y = n, fill = country)) +
  geom_bar(stat = "identity", width = 1, linewidth = 0.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(percent,1),"%")), 
            position = position_stack(vjust = 0.5), 
            colour = "black", size = 4.2) +
  scale_fill_manual(values = country_colors) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 11, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, colour = "black")
  ) 

print(p_country)

######2  sex prop pie
#sex_ration<-c(1,1,2,1,1,2,1,1,1,2,2,2,1,2,1,2,2,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,1)
sex_ration<-metainfo$sex
sex_df<-data.frame(sex=sex_ration)|>count(sex,name="n")|>mutate(group=case_when(sex=="f"~"Female",sex=="m"~"Male"),percent=n/sum(n)*100)
sex_col<-c("Male"="#6699CC","Female"="#86C7B4")
p_sex<-ggplot(sex_df,aes(x="",y=n,fill=group))+geom_bar(stat="identity",width=1,linewidth=0.5)+coord_polar("y",start=0)+geom_text(aes(label=paste0(round(percent,1),"%")),position=position_stack(vjust=0.5),colour="black",size=4.2)+scale_fill_manual(values=sex_col)+theme_void()+theme(legend.title=element_blank(),legend.text=element_text(size=11,colour="black"),plot.title=element_text(hjust=0.5,size=14,colour="black"))
print(p_sex)

table(sex_ration)
#Female-to-male ratio
17/19

######sup  pathway prop pie
pathway_col<-c(
  "NFKB"="#B3D1E7",
  "Inflammasome"="#B07A99",
  "Immune metabolism"="#847AB3",
  "Interferon regulation"="#EC706E",
  "IECD"="#FEC260",
  "Uncategoried"="#EBB1A4")


pathway_df <- metainfo %>% count(pathway, name = "n") %>% mutate(percent = n / sum(n) * 100) %>% arrange(desc(percent))
pathway_df$pathway<-factor(pathway_df$pathway,levels = pathway_df$pathway)
p_pathway<-ggplot(pathway_df,aes(x="",y=n,fill=pathway))+geom_bar(stat="identity",width=1,linewidth=0.5)+coord_polar("y",start=0)+geom_text(aes(label=paste0(round(percent,1),"%")),position=position_stack(vjust=0.5),colour="black",size=4.2)+scale_fill_manual(values=pathway_col)+theme_void()+theme(legend.title=element_blank(),legend.text=element_text(size=11,colour="black"),plot.title=element_text(hjust=0.5,size=14,colour="black"))
table(metainfo$pathway)/nrow(metainfo)

library(patchwork)
p_country|p_sex|p_pathway


######3 Bar plot for mutation type, novel or reported and genotype
library(tidyverse)
library(patchwork)
plot_meta_bar <- function(metainfo, target_col, fill_color ="#6699CC", 
                  bar_width = 0.6) {
  metainfo1<-metainfo[!is.na(metainfo[[target_col]]),]
  split_vec <- metainfo1[[target_col]] %>%
    str_split("/") %>% unlist() %>% str_trim()
  cnt <- tibble(val = split_vec) %>%
    count(val) %>%
    mutate(total = sum(n), percent = n / total * 100) %>%
    arrange(desc(percent))
  cnt$val <- factor(cnt$val, levels = cnt$val)
  p <- ggplot(cnt, aes(x = val, y = percent)) +
    geom_col(fill = fill_color, width = bar_width) +
    geom_text(aes(label = round(percent, 2)), vjust = -0.25, size = 3.5, colour = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
      axis.ticks.x = element_line(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks.y = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(colour = "black")
    ) +
    labs(y = "Percentage (%)")
  return(p)
}

p_mut <- plot_meta_bar(metainfo, "MUT")
p_geno <- plot_meta_bar(metainfo, "genotype") + theme(axis.title.y = element_blank())
p_novel <- plot_meta_bar(metainfo, "novel") + theme(axis.title.y = element_blank())
p_acmg <- plot_meta_bar(metainfo, "ACMG_class") + theme(axis.title.y = element_blank())
p_mecha <- plot_meta_bar(metainfo, "mechanism") + theme(axis.title.y = element_blank())

p_mut | p_geno | p_novel | p_mecha |p_acmg

table(metainfo$MUT)

# Frameshift Frameshift/Missense 
# 4                   1 
# Large deletion            Missense 
# 1                  18 
# Missense/Frameshift   Missense/Missense 
# 3                   2 
# Nonsense Nonsense/Frameshift 
# 5                   1 
# Splicing 
# 1 


# table(split_vec)
# split_vec
# Frameshift Large deletion       Missense 
# 9              1             26 
# Nonsense       Splicing 
# 6              1 

# table(metainfo$novel)
# 
# novel variant novel variant in known site 
# 4                           2 
# reported           reported/reported 
# 24                           6 
# table(  split_vec)
# split_vec
# novel variant novel variant in known site 
# 4                           2 
# reported 
# 36 

##check
MUT1<-c()
for(i in 1:nrow(metainfo)){
  length(grep("\\/",metainfo$MUT[i]))->l
  if(length(l)!=0){
    res<-unlist(strsplit(metainfo$MUT[i],split = "/"))
  }else{
    res<-metainfo$MUT[i]
  }
  MUT1<-c(MUT1,res)
}
table(MUT1)/sum(table(MUT1))

novel1<-c()
for(i in 1:nrow(metainfo)){
  length(grep("\\/",metainfo$novel[i]))->l
  if(length(l)!=0){
    res<-unlist(strsplit(metainfo$novel[i],split = "/"))
  }else{
    res<-metainfo$novel[i]
  }
  novel1<-c(novel1,res)
}
table(novel1)/sum(table(novel1))

ACMG_class1<-c()
for(i in 1:nrow(metainfo)){
  length(grep("\\/",metainfo$ACMG_class[i]))->l
  if(length(l)!=0){
    res<-unlist(strsplit(metainfo$ACMG_class[i],split = "/"))
  }else{
    res<-metainfo$ACMG_class[i]
  }
  ACMG_class1<-c(ACMG_class1,res)
}
table(ACMG_class1)/sum(table(ACMG_class1))


######4 age plot
setwd("E:\\Cohort PPT\\JIA")
library(openxlsx)
clinicalinfo<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
metainfo<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")
####(1) scatter distribution group by different mechanism
# distri_age v2: onset vs diagnosis scatter faceted by mechanism
distri_age <- function(gene_colors, abline_lwd = 0.4, point_size = 3, text_size = 2.8) {
  library(openxlsx); library(ggplot2); library(dplyr); library(patchwork); library(ggrepel)
  gene_order <- c("ADA2","LACC1","ELF4","XIAP","MVK","LPIN2","NLRP3","NOD2",
                  "PSTPIP1","RAC2","COPA","UNC93B1","SOCS1","TNFAIP3")
  if (length(gene_colors) != length(gene_order))
    warning("gene_colors length (", length(gene_colors), ") != gene_order length (", length(gene_order), ")")
  if (is.null(names(gene_colors))) names(gene_colors) <- gene_order
  cat("Color mapping:\n"); print(data.frame(gene=names(gene_colors), color=gene_colors, row.names=NULL))
  
  metainfo <- read.xlsx("0903_Clean case.xlsx")
  metainfo$ID <- trimws(metainfo$ID); metainfo$gene <- toupper(trimws(metainfo$gene))
  metainfo$mechanism <- tolower(trimws(metainfo$mechanism))
  metainfo$mechanism[is.na(metainfo$mechanism)|metainfo$mechanism==""] <- "na"
  metainfo$sex <- tolower(trimws(metainfo$sex))
  metainfo$sex[metainfo$sex=="m"] <- "male"; metainfo$sex[metainfo$sex=="f"] <- "female"
  
  clin <- read.xlsx("code\\Clinical data.xlsx", colNames=FALSE)
  pids <- trimws(as.character(clin[1, 3:ncol(clin)]))
  extr <- function(x) { x <- trimws(as.character(x)); x[grepl("not available|^$|^na$",x,ignore.case=TRUE)] <- NA; as.numeric(gsub("[^0-9.]","",x)) }
  age_df <- data.frame(ID=pids,
                       onset=extr(clin[grep("onset",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]),
                       diag =extr(clin[grep("diagnosis",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]),
                       stringsAsFactors=FALSE)
  df <- merge(metainfo[,c("ID","gene","mechanism","sex")], age_df, by="ID", all.x=TRUE)
  df <- df[df$mechanism != "na" & !is.na(df$onset), ]
  cat("Patients after filtering (excl. TNFRSF1A/na):", nrow(df), "\n")
  
  mech_list <- list("Deficiency"="deficiency","Gain-of-function"="gain-of-function","Haploinsufficiency"="haploinsufficiency")
  plot_list <- list()
  for (i in seq_along(mech_list)) {
    mn <- names(mech_list)[i]; sub <- df[df$mechanism == mech_list[[mn]], ]
    g_in <- unique(sub$gene)
    max_xy <- max(sub$onset, sub$diag, na.rm=TRUE) * 1.1
    p <- ggplot(sub, aes(x=onset, y=diag)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="#661E01", linewidth=abline_lwd) +
      geom_point(aes(color=gene, shape=sex), size=point_size, stroke=0.5) +
      geom_text_repel(aes(label=ID), size=text_size, color="black",
                      segment.color="grey50", segment.size=0.3, segment.alpha=0.7,
                      box.padding=0.3, point.padding=0.2, max.overlaps=50, seed=42) +
      scale_color_manual(values=gene_colors[g_in], breaks=g_in) +
      scale_shape_manual(values=c("male"=16,"female"=17)) +
      coord_cartesian(xlim=c(0, max_xy), ylim=c(0, max_xy)) +
      labs(title=mn, x="Age of onset (y)", y="Age of diagnosis (y)", color=NULL, shape=NULL) +
      theme_classic() +
      theme(
        panel.border      = element_rect(color="black", fill=NA, linewidth=0.5),
        axis.line         = element_blank(),
        panel.grid.major  = element_line(color="grey85", linetype="dashed", linewidth=0.3),
        panel.grid.minor  = element_line(color="grey92", linetype="dashed", linewidth=0.2),
        axis.text         = element_text(color="black", size=9, face="plain"),
        axis.title        = element_text(color="black", size=10, face="plain"),
        plot.title        = element_text(size=11, face="plain", hjust=0.5),
        legend.text       = element_text(size=8, face="plain"),
        legend.key.size   = unit(0.35,"cm"),
        legend.position   = "bottom",
        legend.box        = "vertical",
        legend.margin     = margin(0,0,0,0),
        plot.margin       = if (i==1) margin(8,4,8,8) else margin(8,4,8,0)
      )
    if (i > 1) p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    plot_list[[mn]] <- p
  }
  wrap_plots(plot_list, ncol=3) + plot_layout(guides="collect") & theme(legend.position="bottom")
}
my_colors <- c("#0B71AB","#6699CC","#B3D1E7","#5C6BC0","#A992C0","#C8BBD5",
               "#FFEFC1","#FEC260","#F49D5C","#C95968","#EBB1A4","#8C2522",
               "#96AF95","#008280")
distri_age(gene_colors = my_colors, abline_lwd = 0.5, point_size = 3, text_size = 2.8)
##(a)median and IQR age for three groups
calc_age_stats <- function() {
  library(openxlsx); library(dplyr)
  metainfo <- read.xlsx("0903_Clean case.xlsx")
  metainfo$ID <- trimws(metainfo$ID); metainfo$gene <- toupper(trimws(metainfo$gene))
  metainfo$mechanism <- tolower(trimws(metainfo$mechanism))
  metainfo$mechanism[is.na(metainfo$mechanism)|metainfo$mechanism==""] <- "na"
  clin <- read.xlsx("code\\Clinical data.xlsx", colNames=FALSE)
  pids <- trimws(as.character(clin[1, 3:ncol(clin)]))
  extr <- function(x) { x <- trimws(as.character(x)); x[grepl("not available|^$|^na$",x,ignore.case=TRUE)] <- NA; as.numeric(gsub("[^0-9.]","",x)) }
  age_df <- data.frame(ID=pids, onset=extr(clin[grep("onset",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]), diag=extr(clin[grep("diagnosis",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]), stringsAsFactors=FALSE)
  df <- merge(metainfo[,c("ID","gene","mechanism")], age_df, by="ID", all.x=TRUE)
  df <- df[df$mechanism != "na" & !is.na(df$onset), ]
  mech_list <- c("Deficiency"="deficiency","Gain-of-function"="gain-of-function","Haploinsufficiency"="haploinsufficiency")
  do.call(rbind, lapply(seq_along(mech_list), function(i) {
    sub <- df[df$mechanism == mech_list[[i]], ]
    data.frame(mechanism=names(mech_list)[i], n=nrow(sub),
               onset_median=round(median(sub$onset,na.rm=TRUE),1),
               onset_IQR=paste0(round(quantile(sub$onset,0.25,na.rm=TRUE),1),"-",round(quantile(sub$onset,0.75,na.rm=TRUE),1)),
               diag_median=round(median(sub$diag,na.rm=TRUE),1),
               diag_IQR=paste0(round(quantile(sub$diag,0.25,na.rm=TRUE),1),"-",round(quantile(sub$diag,0.75,na.rm=TRUE),1)),
               stringsAsFactors=FALSE)
  }))
}
age_stats <- calc_age_stats()
print(age_stats)

####(2)forest and IQR age ， foreset plot
plot_age_forest_both <- function(onset_color, diag_color, point_size=3.5, line_width=0.6, vline_width=0.4) {
  library(openxlsx); library(ggplot2); library(dplyr); library(patchwork)
  metainfo <- read.xlsx("0903_Clean case.xlsx")
  metainfo$ID <- trimws(metainfo$ID); metainfo$gene <- toupper(trimws(metainfo$gene))
  metainfo$mechanism <- tolower(trimws(metainfo$mechanism))
  metainfo$mechanism[is.na(metainfo$mechanism)|metainfo$mechanism==""] <- "na"
  clin <- read.xlsx("code\\Clinical data.xlsx", colNames=FALSE)
  pids <- trimws(as.character(clin[1, 3:ncol(clin)]))
  extr <- function(x) { x <- trimws(as.character(x)); x[grepl("not available|^$|^na$",x,ignore.case=TRUE)] <- NA; as.numeric(gsub("[^0-9.]","",x)) }
  age_df <- data.frame(ID=pids, onset=extr(clin[grep("onset",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]), diag=extr(clin[grep("diagnosis",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]), stringsAsFactors=FALSE)
  df <- merge(metainfo[,c("ID","gene","mechanism")], age_df, by="ID", all.x=TRUE)
  df <- df[df$mechanism != "na" & !is.na(df$onset), ]
  mech_list <- c("Deficiency"="deficiency","Gain-of-function"="gain-of-function","Haploinsufficiency"="haploinsufficiency")
  stats_list <- list(); delay_list <- list()
  for (i in seq_along(mech_list)) {
    mn <- names(mech_list)[i]; sub <- df[df$mechanism == mech_list[[i]], ]
    stats_list[[length(stats_list)+1]] <- data.frame(mechanism=mn, type="Onset", median=round(median(sub$onset,na.rm=TRUE),2), q25=round(as.numeric(quantile(sub$onset,0.25,na.rm=TRUE)),2), q75=round(as.numeric(quantile(sub$onset,0.75,na.rm=TRUE)),2), n=nrow(sub))
    stats_list[[length(stats_list)+1]] <- data.frame(mechanism=mn, type="Diagnosis", median=round(median(sub$diag,na.rm=TRUE),2), q25=round(as.numeric(quantile(sub$diag,0.25,na.rm=TRUE)),2), q75=round(as.numeric(quantile(sub$diag,0.75,na.rm=TRUE)),2), n=nrow(sub))
    delay_list[[length(delay_list)+1]] <- data.frame(mechanism=mn, n=nrow(sub), onset_median=round(median(sub$onset,na.rm=TRUE),2), diag_median=round(median(sub$diag,na.rm=TRUE),2), delay_median=round(median(sub$diag,na.rm=TRUE)-median(sub$onset,na.rm=TRUE),2))
  }
  stats_df <- do.call(rbind, stats_list)
  delay_df <- do.call(rbind, delay_list)
  stats_df$type <- factor(stats_df$type, levels=c("Onset","Diagnosis"))
  stats_df$mechanism <- factor(stats_df$mechanism, levels=c("Deficiency","Gain-of-function","Haploinsufficiency"))
  base <- ggplot(stats_df, aes(x=median, y=type, color=type)) +
    geom_vline(aes(xintercept=median), linetype="dashed", color="gray90", linewidth=vline_width) +
    geom_errorbarh(aes(xmin=q25, xmax=q75, height=0.25), linewidth=line_width) +
    geom_point(size=point_size) +
    geom_text(aes(label=paste0(sprintf("%.2f", median), " (n=", n, ")")), hjust=-0.3, vjust=0.5, size=2.8, color="black") +
    scale_color_manual(values=c("Onset"=onset_color, "Diagnosis"=diag_color)) +
    scale_x_continuous(expand=expansion(mult=c(0.1,0.25))) +
    labs(x="Age (years)", y=NULL, color=NULL) +
    theme_classic() +
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.5),
      axis.line = element_blank(),
      axis.text = element_text(color="black", size=9),
      axis.title = element_text(color="black", size=10),
      strip.text = element_text(size=10, face="plain"),
      strip.background = element_rect(color="black", linewidth=0.5),
      legend.position = "bottom",
      legend.key.size = unit(0.35,"cm"),
      plot.margin = margin(8,8,8,8)
    )
  p_v <- base + facet_grid(mechanism ~ ., scales="free_y", space="free_y")
  p_h <- base + facet_grid(. ~ mechanism, scales="free_x")
  list(vertical=p_v, horizontal=p_h, stats=stats_df, delay=delay_df)
}

res <- plot_age_forest_both(onset_color="#6699CC", diag_color="#D8A0C0")
res$vertical
res$horizontal
res$stats
res$delay

##(a) check the proportion of patients with diagnosed age larger than onset
calc_delay_ratio <- function() {
  library(openxlsx); library(dplyr)
  metainfo <- read.xlsx("0903_Clean case.xlsx")
  metainfo$ID <- trimws(metainfo$ID); metainfo$gene <- toupper(trimws(metainfo$gene))
  metainfo$mechanism <- tolower(trimws(metainfo$mechanism))
  metainfo$mechanism[is.na(metainfo$mechanism)|metainfo$mechanism==""] <- "na"
  clin <- read.xlsx("code\\Clinical data.xlsx", colNames=FALSE)
  pids <- trimws(as.character(clin[1, 3:ncol(clin)]))
  extr <- function(x) { x <- trimws(as.character(x)); x[grepl("not available|^$|^na$",x,ignore.case=TRUE)] <- NA; as.numeric(gsub("[^0-9.]","",x)) }
  age_df <- data.frame(ID=pids, onset=extr(clin[grep("onset",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]), diag=extr(clin[grep("diagnosis",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]), stringsAsFactors=FALSE)
  df <- merge(metainfo[,c("ID","gene","mechanism")], age_df, by="ID", all.x=TRUE)
  df <- df[df$mechanism != "na" & !is.na(df$onset) & !is.na(df$diag), ]
  mech_list <- c("Deficiency"="deficiency","Gain-of-function"="gain-of-function","Haploinsufficiency"="haploinsufficiency")
  do.call(rbind, lapply(seq_along(mech_list), function(i) {
    mn <- names(mech_list)[i]; sub <- df[df$mechanism == mech_list[[i]], ]
    n_total <- nrow(sub)
    n_delayed <- sum(sub$diag > sub$onset, na.rm=TRUE)
    data.frame(mechanism=mn, n_total=n_total, n_diag_gt_onset=n_delayed, pct_diag_gt_onset=round(n_delayed/n_total*100,2), stringsAsFactors=FALSE)
  }))
}

delay_ratio <- calc_delay_ratio()
print(delay_ratio)

#####(3) density plots
distri_age_sex <- function(male_col = "#0072B2", female_col = "#AF478A",
                           med_male_col = "#003D81", med_female_col = "#8C2522",
                           med_overall_col = "#333333", med_lwd = 0.7,
                           line_w = 0.8, fill_alpha = 0.25, text_size = 2.8) {
  library(openxlsx); library(ggplot2); library(dplyr); library(patchwork)
  
  metainfo <- read.xlsx("0903_Clean case.xlsx")
  metainfo$ID <- trimws(metainfo$ID)
  metainfo$sex <- tolower(trimws(metainfo$sex))
  metainfo$sex[metainfo$sex=="m"] <- "male"; metainfo$sex[metainfo$sex=="f"] <- "female"
  
  clin <- read.xlsx("code\\Clinical data.xlsx", colNames=FALSE)
  pids <- trimws(as.character(clin[1, 3:ncol(clin)]))
  extr <- function(x) { x <- trimws(as.character(x)); x[grepl("not available|^$|^na$",x,ignore.case=TRUE)] <- NA; as.numeric(gsub("[^0-9.]","",x)) }
  age_df <- data.frame(ID=pids,
                       onset=extr(clin[grep("onset",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]),
                       diag =extr(clin[grep("diagnosis",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]),
                       stringsAsFactors=FALSE)
  df <- merge(metainfo[,c("ID","sex")], age_df, by="ID", all.x=TRUE)
  df$delay <- df$diag - df$onset
  df <- df[!is.na(df$sex) & df$sex %in% c("male","female") & !is.na(df$onset), ]
  df$sex <- factor(df$sex, levels=c("male","female"))
  sex_colors <- c("male"=male_col, "female"=female_col)
  
  kw <- list(onset=kruskal.test(onset~sex,df), diag=kruskal.test(diag~sex,df), delay=kruskal.test(delay~sex,df))
  cat("Kruskal-Wallis p -- onset:",round(kw$onset$p.value,4)," diag:",round(kw$diag$p.value,4)," delay:",round(kw$delay$p.value,4),"\n")
  
  med <- c(male=median(df$onset[df$sex=="male"],na.rm=TRUE),
           female=median(df$onset[df$sex=="female"],na.rm=TRUE),
           overall=median(df$onset,na.rm=TRUE))
  cat("Onset median -- male:",med["male"]," female:",med["female"]," overall:",med["overall"],"\n")
  
  p_lab <- function(pv) ifelse(pv<0.001,"p < 0.001",paste0("p = ",round(pv,3)))
  
  base_th <- theme_classic() + theme(
    panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
    axis.line=element_blank(),
    panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    axis.text=element_text(color="black",size=9,face="plain"),
    axis.title=element_text(color="black",size=10,face="plain"),
    plot.title=element_text(size=11,face="plain",hjust=0.5),
    legend.text=element_text(size=8,face="plain"), legend.key.size=unit(0.35,"cm"),
    legend.position="bottom")
  
  # --- Plot 1: onset, median lines + top-right legend, KW p top-left ---
  p1 <- ggplot(df, aes(x=onset, fill=sex, color=sex)) +
    geom_density(linewidth=line_w, alpha=fill_alpha) +
    geom_vline(xintercept=med["male"],    color=med_male_col,    linetype="dashed", linewidth=med_lwd) +
    geom_vline(xintercept=med["female"],  color=med_female_col,  linetype="dashed", linewidth=med_lwd) +
    geom_vline(xintercept=med["overall"], color=med_overall_col, linetype="dotted", linewidth=med_lwd) +
    annotate("text", x=Inf, y=Inf, label=paste0("Male median: ",   med["male"],    " y"),
             color=med_male_col,    size=text_size*0.85, hjust=1.1, vjust=1.6) +
    annotate("text", x=Inf, y=Inf, label=paste0("Female median: ", med["female"],  " y"),
             color=med_female_col,  size=text_size*0.85, hjust=1.1, vjust=3.0) +
    annotate("text", x=Inf, y=Inf, label=paste0("Overall median: ",med["overall"], " y"),
             color=med_overall_col, size=text_size*0.85, hjust=1.1, vjust=4.4) +
    annotate("text", x=-Inf, y=Inf, label=paste0("Kruskal-Wallis\n",p_lab(kw$onset$p.value)),
             hjust=-0.1, vjust=1.3, size=text_size, color="black") +
    scale_fill_manual(values=sex_colors, labels=c("Male","Female")) +
    scale_color_manual(values=sex_colors, labels=c("Male","Female")) +
    labs(title="Age of onset", x="Years", y="Density", fill=NULL, color=NULL) +
    base_th + theme(plot.margin=margin(8,4,8,8))
  
  # --- Plot 2: diagnosis age ---
  p2 <- ggplot(df, aes(x=diag, fill=sex, color=sex)) +
    geom_density(linewidth=line_w, alpha=fill_alpha) +
    annotate("text", x=Inf, y=Inf, label=paste0("Kruskal-Wallis\n",p_lab(kw$diag$p.value)), hjust=1.1, vjust=1.3, size=text_size, color="black") +
    scale_fill_manual(values=sex_colors) + scale_color_manual(values=sex_colors) +
    labs(title="Age of diagnosis", x="Years", y=NULL, fill=NULL, color=NULL) +
    base_th + theme(plot.margin=margin(8,4,8,0), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  # --- Plot 3: diagnostic delay ---
  p3 <- ggplot(df, aes(x=delay, fill=sex, color=sex)) +
    geom_density(linewidth=line_w, alpha=fill_alpha) +
    annotate("text", x=Inf, y=Inf, label=paste0("Kruskal-Wallis\n",p_lab(kw$delay$p.value)), hjust=1.1, vjust=1.3, size=text_size, color="black") +
    scale_fill_manual(values=sex_colors) + scale_color_manual(values=sex_colors) +
    labs(title="Diagnostic delay", x="Years", y=NULL, fill=NULL, color=NULL) +
    base_th + theme(plot.margin=margin(8,4,8,0), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  wrap_plots(p1, p2, p3, ncol=3) + plot_layout(guides="collect") & theme(legend.position="bottom")
}

age_density3plots<-distri_age_sex(
  male_col="#0072B2", female_col="#AF478A",
  med_male_col="#003D81", med_female_col="#8C2522", med_overall_col="#333333",
  med_lwd=0.7, line_w=0.8, fill_alpha=0.25, text_size=2.8)
age_density3plots

##(a) check median onset age
Clinical_data<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
age_data<-Clinical_data[
Clinical_data$Number%in%c("Sex","Age of onset (year)","Age of diagnosis (year)"),]
rownames(age_data)<-age_data$Number
age_data<-age_data[,-c(1,2)]
age_data<-as.data.frame(t(age_data))
age_data$`Age of onset (year)`<-as.numeric(gsub("y","",age_data$`Age of onset (year)`))
age_data$`Age of diagnosis (year)`<-as.numeric(gsub("y","",age_data$`Age of diagnosis (year)`))
age_data$`Age of delay (year)`<-age_data$`Age of diagnosis (year)`-age_data$`Age of onset (year)`

median(age_data[age_data$Sex=="female","Age of onset (year)"])#2
median(age_data[age_data$Sex=="male","Age of onset (year)"])#2.085
median(age_data$`Age of onset (year)`)#2 

##(b) check  three wilcoxon test P value
age_check <- age_data[!is.na(age_data$`Age of onset (year)`), ]
age_check$Sex <- tolower(trimws(age_check$Sex))
kw_check <- list(
  onset = kruskal.test(`Age of onset (year)` ~ Sex, data = age_check),
  diag  = kruskal.test(`Age of diagnosis (year)` ~ Sex, data = age_check),
  delay = kruskal.test(`Age of delay (year)` ~ Sex, data = age_check)
)
cat("Kruskal-Wallis p -- onset:", round(kw_check$onset$p.value,4),
    " diag:", round(kw_check$diag$p.value,4),
    " delay:", round(kw_check$delay$p.value,4), "\n")


####(4) scatter with age value and median, IQR,  group with genes and mechanism
metainfo <- read.xlsx("0903_Clean case.xlsx")
metainfo$ID <- trimws(metainfo$ID)
metainfo$gene <- toupper(trimws(metainfo$gene))
metainfo$mechanism <- tolower(trimws(metainfo$mechanism))
metainfo$mechanism[is.na(metainfo$mechanism)|metainfo$mechanism==""] <- "na"

clin <- read.xlsx("code\\Clinical data.xlsx", colNames=FALSE)
pids <- trimws(as.character(clin[1, 3:ncol(clin)]))
extr <- function(x) { x <- trimws(as.character(x)); x[grepl("not available|^$|^na$",x,ignore.case=TRUE)] <- NA; as.numeric(gsub("[^0-9.]","",x)) }
age_df <- data.frame(ID=pids,
                     onset=extr(clin[grep("onset",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]),
                     diag =extr(clin[grep("diagnosis",clin[,1],ignore.case=TRUE)[1], 3:ncol(clin)]),
                     stringsAsFactors=FALSE)
age_df$delay <- age_df$diag - age_df$onset

df <- merge(metainfo[,c("ID","gene","mechanism")], age_df, by="ID", all.x=TRUE)
df <- df[!is.na(df$onset), ]
cat("Patients:", nrow(df), "\n")
print(sort(table(df$gene)))

plot_gene_age <- function(data, age_col, ylab, ybreaks,
                          point_color = "#6699CC", bar_color = "#BB0021",
                          point_size = 1.5, bar_width = 0.4, bar_lwd = 0.5,
                          mech_colors = c("deficiency"="#6699CC",
                                          "gain-of-function"="#FEC260",
                                          "haploinsufficiency"="#AF478A",
                                          "na"="black")) {
  d <- data[!is.na(data[[age_col]]), ]
  norm_p <- shapiro.test(d[[age_col]])$p.value
  kw_p <- kruskal.test(d[[age_col]] ~ d$gene)$p.value
  cat("Normality P =", round(norm_p,4), "| Kruskal-Wallis P =", round(kw_p,4), "\n")
  gs <- d %>% group_by(gene) %>%
    summarise(median_age=median(.data[[age_col]]),
              q1=quantile(.data[[age_col]],0.25),
              q3=quantile(.data[[age_col]],0.75), .groups="drop") %>%
    arrange(median_age)
  gs$mech <- sapply(gs$gene, function(g) unique(d$mechanism[d$gene==g])[1])
  gs$lab_col <- mech_colors[gs$mech]
  gs$lab_col[is.na(gs$lab_col)] <- "black"
  print(gs)
  d$gene <- factor(d$gene, levels=gs$gene)
  gs$gene <- factor(gs$gene, levels=gs$gene)
  ggplot(d, aes(x=gene, y=.data[[age_col]])) +
    geom_point(color=point_color, size=point_size) +
    geom_errorbar(data=gs, aes(y=median_age, ymin=q1, ymax=q3), color=bar_color, width=bar_width, linewidth=bar_lwd) +
    geom_point(data=gs, aes(y=median_age), color=bar_color, size=2, shape=16) +
    scale_y_continuous(breaks=ybreaks, expand=expansion(mult=c(0.05,0.1))) +
    theme_classic() +
    xlab(NULL) + ylab(ylab) +
    theme(axis.text=element_text(color="black", size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=10, colour=gs$lab_col),
          axis.title=element_text(color="black"))
}

plot_gene_age_mech <- function(data, age_col, ylab, ybreaks,
                               point_color = "#6699CC", bar_color = "#BB0021",
                               point_size = 1.5, bar_width = 0.4, bar_lwd = 0.5,
                               mech_order = c("na","gain-of-function","haploinsufficiency","deficiency")) {
  d <- data[!is.na(data[[age_col]]), ]
  norm_p <- shapiro.test(d[[age_col]])$p.value
  kw_p <- kruskal.test(d[[age_col]] ~ d$gene)$p.value
  cat("Normality P =", round(norm_p,4), "| Kruskal-Wallis P =", round(kw_p,4), "\n")
  gs <- d %>% group_by(gene) %>%
    summarise(n=n(), median_age=median(.data[[age_col]]),
              q1=quantile(.data[[age_col]],0.25),
              q3=quantile(.data[[age_col]],0.75), .groups="drop")
  gs$mech <- sapply(gs$gene, function(g) unique(d$mechanism[d$gene==g])[1])
  gs$mech[is.na(gs$mech)] <- "na"
  gs$mech <- factor(gs$mech, levels=mech_order)
  gs <- gs %>% arrange(mech, median_age)
  print(gs)
  d$gene <- factor(d$gene, levels=gs$gene)
  gs$gene <- factor(gs$gene, levels=gs$gene)
  ge <- gs[gs$n >= 2, ]
  p <- ggplot(d, aes(x=gene, y=.data[[age_col]])) +
    geom_point(color=point_color, size=point_size) +
    scale_y_continuous(breaks=ybreaks, expand=expansion(mult=c(0.05,0.1))) +
    theme_classic() +
    xlab(NULL) + ylab(ylab) +
    theme(axis.text=element_text(color="black", size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.title=element_text(color="black"))
  if (nrow(ge) > 0)
    p <- p +
    geom_errorbar(data=ge, aes(y=median_age, ymin=q1, ymax=q3), color=bar_color, width=bar_width, linewidth=bar_lwd) +
    geom_point(data=ge, aes(y=median_age), color=bar_color, size=2, shape=16)
  p
}

plot_gene_age_combo <- function(data, age_col, ylab, ybreaks,
                                point_color = "#6699CC", bar_color = "#BB0021",
                                point_size = 1.5, bar_width = 0.4, bar_lwd = 0.5,
                                mech_colors = c("deficiency"="#6699CC",
                                                "gain-of-function"="#FEC260",
                                                "haploinsufficiency"="#AF478A",
                                                "na"="black"),
                                mech_order = c("gain-of-function","haploinsufficiency","deficiency"),
                                left_width = 0.9, right_width = 0.4,
                                jitter_w = 0.15, seed = 42) {
  d <- data[!is.na(data[[age_col]]), ]
  norm_p <- shapiro.test(d[[age_col]])$p.value
  kw_p <- kruskal.test(d[[age_col]] ~ d$gene)$p.value
  cat("Normality P =", round(norm_p,4), "| Kruskal-Wallis P =", round(kw_p,4), "\n")
  gs <- d %>% group_by(gene) %>%
    summarise(n=n(), median_age=median(.data[[age_col]]),
              q1=quantile(.data[[age_col]],0.25),
              q3=quantile(.data[[age_col]],0.75), .groups="drop")
  gs$mech <- sapply(gs$gene, function(g) unique(d$mechanism[d$gene==g])[1])
  gs$mech[is.na(gs$mech)] <- "na"
  gs$mech <- factor(gs$mech, levels=c("na", mech_order))
  gs <- gs %>% arrange(mech, median_age)
  print(gs)
  d$gene <- factor(d$gene, levels=gs$gene)
  gs$gene <- factor(gs$gene, levels=gs$gene)
  ge <- gs[gs$n >= 2, ]
  p_left <- ggplot(d, aes(x=gene, y=.data[[age_col]])) +
    geom_point(color=point_color, size=point_size) +
    scale_y_continuous(breaks=ybreaks, expand=expansion(mult=c(0.05,0.1))) +
    theme_classic() +
    xlab(NULL) + ylab(ylab) +
    theme(axis.text=element_text(color="black", size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.title=element_text(color="black"))
  if (nrow(ge) > 0)
    p_left <- p_left +
    geom_errorbar(data=ge, aes(y=median_age, ymin=q1, ymax=q3), color=bar_color, width=bar_width, linewidth=bar_lwd) +
    geom_point(data=ge, aes(y=median_age), color=bar_color, size=2, shape=16)
  dm <- d[d$mechanism %in% mech_order, ]
  dm$mechanism <- factor(dm$mechanism, levels=mech_order)
  ms <- dm %>% group_by(mechanism) %>%
    summarise(n=n(), median_age=median(.data[[age_col]]),
              q1=quantile(.data[[age_col]],0.25),
              q3=quantile(.data[[age_col]],0.75), .groups="drop")
  print(ms)
  p_right <- ggplot(dm, aes(x=mechanism, y=.data[[age_col]])) +
    geom_point(aes(color=mechanism), size=point_size,
               position=position_jitter(width=jitter_w, height=0, seed=seed)) +
    geom_errorbar(data=ms, aes(y=median_age, ymin=q1, ymax=q3), color=bar_color, width=bar_width*0.5, linewidth=bar_lwd) +
    geom_point(data=ms, aes(y=median_age), color=bar_color, size=2, shape=16) +
    scale_color_manual(values=mech_colors[mech_order]) +
    scale_y_continuous(breaks=ybreaks, expand=expansion(mult=c(0.05,0.1))) +
    theme_classic() +
    xlab(NULL) + ylab(NULL) +
    theme(axis.text=element_text(color="black", size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=10),
          axis.title=element_text(color="black"),
          legend.position="none")
  p_left + p_right + plot_layout(widths=c(left_width, right_width))
}

age_onset_plot  <- plot_gene_age(df, "onset", "Age at onset", c(0,2,6,10,15,20))
age_delay_plot  <- plot_gene_age(df, "delay", "Diagnostic delay", seq(0,35,5))
age_onset_plot2 <- plot_gene_age_mech(df, "onset", "Age at onset", c(0,2,6,10,15,20))
age_delay_plot2 <- plot_gene_age_mech(df, "delay", "Diagnostic delay", seq(0,35,5))
age_onset_combo <- plot_gene_age_combo(df, "onset", "Age at onset", c(0,2,6,10,15,20))
age_delay_combo <- plot_gene_age_combo(df, "delay", "Diagnostic delay", seq(0,35,5))
age_onset_plot
age_delay_plot
age_onset_plot2
age_delay_plot2
age_onset_combo/age_delay_combo


##(a) check and output for manuscript
Clinical_data<-read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx")
age_data<-Clinical_data[
  Clinical_data$Number%in%c("Sex","Age of onset (year)","Age of diagnosis (year)"),]
rownames(age_data)<-age_data$Number
age_data<-age_data[,-c(1,2)]
age_data<-as.data.frame(t(age_data))
age_data$`Age of onset (year)`<-as.numeric(gsub("y","",age_data$`Age of onset (year)`))
age_data$`Age of diagnosis (year)`<-as.numeric(gsub("y","",age_data$`Age of diagnosis (year)`))
age_data$`Age of delay (year)`<-age_data$`Age of diagnosis (year)`-age_data$`Age of onset (year)`

age_data$mechanism<-metainfo[match(rownames(age_data),metainfo$ID),'mechanism']
age_data$gene<-metainfo[match(rownames(age_data),metainfo$ID),'gene']


summarize_age_stats <- function(data, age_col,
                                mech_order = c("gain-of-function","haploinsufficiency","deficiency")) {
  library(dplyr)
  d <- data[!is.na(data[[age_col]]) & !is.na(data$gene), ]
  gs <- d %>% group_by(gene) %>%
    summarise(n=n(),
              median=median(.data[[age_col]]),
              q1=quantile(.data[[age_col]],0.25),
              q3=quantile(.data[[age_col]],0.75),
              .groups="drop") %>%
    mutate(iqr=q3-q1)
  gs[gs$n==1, c("q1","q3","iqr")] <- NA
  gs <- gs %>% arrange(desc(n), median)
  dm <- d[!is.na(d$mechanism) & trimws(d$mechanism)!="" & tolower(trimws(d$mechanism))!="na", ]
  dm$mechanism <- tolower(trimws(dm$mechanism))
  ms <- dm %>% group_by(mechanism) %>%
    summarise(n=n(),
              median=median(.data[[age_col]]),
              q1=quantile(.data[[age_col]],0.25),
              q3=quantile(.data[[age_col]],0.75),
              .groups="drop") %>%
    mutate(iqr=q3-q1)
  lev <- intersect(c(mech_order, "na"), unique(ms$mechanism))
  lev <- c(lev, setdiff(unique(ms$mechanism), lev))
  ms$mechanism <- factor(ms$mechanism, levels=lev)
  ms <- ms %>% arrange(mechanism)
  comps <- combn(unique(as.character(ms$mechanism)), 2, simplify=FALSE)
  wt <- do.call(rbind, lapply(comps, function(pr){
    a <- dm[[age_col]][dm$mechanism==pr[1]]
    b <- dm[[age_col]][dm$mechanism==pr[2]]
    pw <- tryCatch(wilcox.test(a,b)$p.value, error=function(e) NA)
    data.frame(group1=pr[1], group2=pr[2], n1=length(a), n2=length(b), p_wilcox=pw)
  }))
  wt$p_signif <- ifelse(is.na(wt$p_wilcox), "ns",
                        ifelse(wt$p_wilcox<0.001,"***",
                               ifelse(wt$p_wilcox<0.01,"**",
                                      ifelse(wt$p_wilcox<0.05,"*","ns"))))
  cat("\n-- by gene --\n"); print(as.data.frame(gs), digits=7)
  cat("\n-- by mechanism --\n"); print(as.data.frame(ms), digits=7)
  cat("\n-- pairwise wilcoxon --\n"); print(as.data.frame(wt), digits=7)
  return(list(gene_stats=gs, mech_stats=ms, wilcox=wt))
}


onset_stats  <- summarize_age_stats(age_data, "Age of onset (year)")
##check again
median(age_data[age_data$mechanism=="gain-of-function"&!is.na(age_data$mechanism),'Age of onset (year)'])
median(age_data[age_data$mechanism=="haploinsufficiency"&!is.na(age_data$mechanism),'Age of onset (year)'])
median(age_data[age_data$mechanism=="deficiency"&!is.na(age_data$mechanism),'Age of onset (year)'])
quantile(age_data[age_data$mechanism=="gain-of-function"&!is.na(age_data$mechanism),'Age of onset (year)'])
quantile(age_data[age_data$mechanism=="haploinsufficiency"&!is.na(age_data$mechanism),'Age of onset (year)'])
quantile(age_data[age_data$mechanism=="deficiency"&!is.na(age_data$mechanism),'Age of onset (year)'])

delay_stats  <- summarize_age_stats(age_data, "Age of delay (year)")
##check again
median(age_data[age_data$mechanism=="gain-of-function"&!is.na(age_data$mechanism),'Age of delay (year)'])
median(age_data[age_data$mechanism=="haploinsufficiency"&!is.na(age_data$mechanism),'Age of delay (year)'])
median(age_data[age_data$mechanism=="deficiency"&!is.na(age_data$mechanism),'Age of delay (year)'])
quantile(age_data[age_data$mechanism=="gain-of-function"&!is.na(age_data$mechanism),'Age of delay (year)'])
quantile(age_data[age_data$mechanism=="haploinsufficiency"&!is.na(age_data$mechanism),'Age of delay (year)'])
quantile(age_data[age_data$mechanism=="deficiency"&!is.na(age_data$mechanism),'Age of delay (year)'])


write_age_stats_xlsx <- function(onset_res, delay_res, file = "age_stats_summary.xlsx"){
  library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, "onset")
  addWorksheet(wb, "delay")
  write_block <- function(sheet, label, df, startRow){
    writeData(wb, sheet, label, startRow = startRow)
    writeData(wb, sheet, df, startRow = startRow + 1)
    startRow + 1 + nrow(df) + 1
  }
  r <- 1
  r <- write_block("onset", "by gene", as.data.frame(onset_res$gene_stats), r)
  r <- write_block("onset", "by mechanism", as.data.frame(onset_res$mech_stats), r)
  write_block("onset", "pairwise wilcoxon", as.data.frame(onset_res$wilcox), r)
  r <- 1
  r <- write_block("delay", "by gene", as.data.frame(delay_res$gene_stats), r)
  r <- write_block("delay", "by mechanism", as.data.frame(delay_res$mech_stats), r)
  write_block("delay", "pairwise wilcoxon", as.data.frame(delay_res$wilcox), r)
  saveWorkbook(wb, file, overwrite = TRUE)
  cat("saved:", file, "\n")
}

write_age_stats_xlsx(onset_stats, delay_stats)

###### gene count for all mutations from previous cases and this study
##(1) filter with domain
library(openxlsx); library(ggplot2); library(dplyr); library(patchwork)
genemeta <- read.xlsx("C:/Users/Administrator/Desktop/15个基因已发表突变中是JIA 表型的搜集/score_meta.xlsx")
f <- c("Autosomal recessive" = "#6699CC",
                  "autosomal dominant" = "#FEC260",
                  "X-linked recessive" = "#AF478A")

plot_gene_bar <- function(genemeta,
                          inher_colors = c(f)){
  library(ggplot2); library(dplyr)
  genemeta <- genemeta[!genemeta$Domain %in% c("\\","NA","") & !is.na(genemeta$Domain), ]
  gene_n <- genemeta %>% count(Gene, name = "n")
  inherit_df <- genemeta %>% group_by(Gene) %>%
    summarise(Inheriance = unique(Inheriance)[1], .groups = "drop") %>%
    left_join(gene_n, by = "Gene")
  inherit_df$Inheriance <- factor(inherit_df$Inheriance,
                                  levels = c("autosomal dominant","Autosomal recessive","X-linked recessive"))
  inherit_df <- inherit_df %>% arrange(Inheriance, desc(n))
  inherit_df$Gene <- factor(inherit_df$Gene, levels = inherit_df$Gene)
  inherit_df$pct <- inherit_df$n / sum(inherit_df$n)
  inherit_df$label <- paste0(round(inherit_df$pct*100,1), "%")
  p <- ggplot(inherit_df, aes(x = Gene, y = pct)) +
    geom_col(aes(fill = Inheriance), width = 0.48, color = NA) +
    geom_text(aes(label = label), vjust = -0.5, size = 3, color = "black") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = inher_colors) +
    facet_grid(~ Inheriance, scales = "free_x", space = "free_x") +
    labs(x = "", y = "Percentage") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          strip.text = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "none")
  return(p)
}
p <- plot_gene_bar(genemeta)
p

##(2) no filter, using this
plot_gene_bar <- function(genemeta,
                          inher_colors = c("Autosomal recessive" = "#6699CC",
                                           "autosomal dominant"   = "#FEC260",
                                           "X-linked recessive"   = "#AF478A")) {
  library(ggplot2); library(dplyr)
  gene_n <- genemeta %>% count(Gene, name = "n")
  inherit_df <- genemeta %>%
    group_by(Gene) %>% summarise(Inheriance = unique(Inheriance)[1], .groups = "drop") %>%
    left_join(gene_n, by = "Gene")
  inherit_df$Inheriance <- factor(inherit_df$Inheriance,
                                  levels = c("autosomal dominant", "Autosomal recessive", "X-linked recessive"))
  inherit_df <- inherit_df %>% arrange(Inheriance, desc(n))
  inherit_df$Gene <- factor(inherit_df$Gene, levels = inherit_df$Gene)
  inherit_df$pct <- inherit_df$n / sum(inherit_df$n)
  inherit_df$label <- paste0(round(inherit_df$pct * 100, 2), "%")
  p <- ggplot(inherit_df, aes(x = Gene, y = pct)) +
    geom_col(aes(fill = Inheriance), width = 0.48, color = NA) +
    geom_text(aes(label = label), vjust = -0.5, size = 3, color = "black") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = inher_colors) +
    facet_grid(~ Inheriance, scales = "free_x", space = "free_x") +
    labs(x = "", y = "Percentage") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          strip.text = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "none")
  return(p)
}
p <- plot_gene_bar(genemeta)
p
##check
prop <- table(genemeta$Gene) / sum(table(genemeta$Gene))
sort(round(prop, 4), decreasing = TRUE) 
sum(round(prop, 4))                     
sum(prop)                                


######6 gene domain for all mutations. Donut plots
####(1) Donut for each gene
library(openxlsx); library(ggplot2); library(dplyr); library(patchwork)
genemeta <- read.xlsx("C:/Users/Administrator/Desktop/15个基因已发表突变中是JIA 表型的搜集/score_meta.xlsx")
genemeta <- genemeta[!genemeta$Domain %in% c("\\","NA","") & !is.na(genemeta$Domain), ]
genemeta$Domain_short <- gsub(" domain", "", genemeta$Domain, ignore.case = TRUE)
gene_n <- genemeta %>% count(Gene, name = "n")
inherit_df <- genemeta %>% group_by(Gene) %>%
  summarise(Inheriance = unique(Inheriance)[1], .groups = "drop") %>%
  left_join(gene_n, by = "Gene")
inherit_df$Inheriance <- factor(inherit_df$Inheriance,
                                levels = c("autosomal dominant","Autosomal recessive","X-linked recessive"))
gene_order <- inherit_df %>% arrange(Inheriance, desc(n)) %>% pull(Gene)
inher_colors <- c("Autosomal recessive" = "#6699CC",
                  "autosomal dominant" = "#FEC260",
                  "X-linked recessive" = "#AF478A")
domain_df <- genemeta %>%
  count(Gene, Domain_short, name = "n") %>%
  group_by(Gene) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  mutate(prop = n / sum(n),
         label = paste0(Domain_short, "\n", round(prop*100,1), "%")) %>%
  ungroup()
domain_palette <- c("#6699CC","#007ABA","#9CD2ED","#86C7B4","#847AB3","#F1BAAB")
plot_list <- list()
for(g in gene_order){
  sub <- domain_df %>% filter(Gene == g) %>% arrange(desc(n))
  n_dom <- nrow(sub)
  sub_colors <- setNames(domain_palette[1:n_dom], sub$Domain_short)
  inh <- inherit_df$Inheriance[inherit_df$Gene == g]
  title_color <- inher_colors[inh]
  if(is.na(title_color)) title_color <- "black"
  p <- ggplot(sub, aes(x = 2, y = prop, fill = Domain_short)) +
    geom_col(width = 1, color = NA) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5),
              size = 2.5, color = "black") +
    scale_fill_manual(values = sub_colors) +
    ggtitle(g) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", color = title_color),
          legend.position = "none")
  plot_list[[g]] <- p
}
p <- wrap_plots(plot_list, nrow = 3, ncol = 5)
p

####(2) Merge all domains of all genes at one bubble plot
library(ggplot2); library(dplyr);library(openxlsx)
genemeta<-read.xlsx("C:\\Users\\Administrator\\Desktop\\15个基因已发表突变中是JIA 表型的搜集\\score_meta.xlsx")
AR  <- c("ADA2","LACC1","MVK","LPIN2")
AD  <- c("TNFRSF1A","NLRP3","COPA","NOD2","PSTPIP1","UNC93B1","TNFAIP3","RAC2","SOCS1")
XLR <- c("ELF4","XIAP")
mech_colors <- c("AR"="#6699CC", "AD"="#FEC260", "XLR"="#AF478A")
mech_order <- c("AR","AD","XLR")


plot_gene_domain_dot <- function(genemeta,
                                 mech_colors = c("AR"="#6699CC","AD"="#FEC260","XLR"="#AF478A"),
                                 mech_order = c("AR","AD","XLR"),
                                 square = FALSE, show_pct = FALSE){
  gm <- genemeta[!is.na(genemeta$Domain) & !genemeta$Domain %in% c("\\","NA",""), ]
  gm$Domain_short <- gsub(" domain", "", gm$Domain, ignore.case = TRUE)
  mech_of <- function(g) factor(ifelse(g %in% AR, "AR", ifelse(g %in% AD, "AD", "XLR")), levels = mech_order)
  gn <- gm %>% count(Gene, name = "n"); gn$mech <- mech_of(gn$Gene)
  gene_order <- gn %>% arrange(mech, desc(n)) %>% pull(Gene)
  df <- gm %>% count(Gene, Domain_short, name = "n") %>% group_by(Gene) %>% mutate(prop = n/sum(n)) %>% ungroup()
  df$Gene <- factor(df$Gene, levels = gene_order)
  df <- df %>% arrange(Gene, desc(n))
  df$dom_id <- factor(paste0(as.character(df$Gene), "_", df$Domain_short),
                      levels = unique(paste0(as.character(df$Gene), "_", df$Domain_short)))
  df$mech <- mech_of(df$Gene)
  dom_lab <- setNames(as.character(df$Domain_short), as.character(df$dom_id))
  block_bounds <- df %>% group_by(Gene) %>% summarise(xmax = max(as.numeric(dom_id)) + 0.5, .groups = "drop")
  mech_sep <- gn %>% group_by(mech) %>% summarise(ymax = max(as.numeric(factor(Gene, levels = gene_order))) + 0.5, .groups = "drop") %>% pull(ymax)
  mech_sep <- mech_sep[-length(mech_sep)]
  p <- ggplot(df, aes(x = dom_id, y = Gene))
  if(square){
    p <- p + geom_point(aes(size = prop, fill = mech), shape = 22, color = "white", linewidth = 0.3, alpha = 0.9)
  } else {
    p <- p + geom_point(aes(size = prop, color = mech), alpha = 0.85)
  }
  if(show_pct)
    p <- p + geom_text(aes(label = paste0(round(prop*100, 1), "%")), size = 2.5, color = "black")
  p <- p +
    scale_x_discrete(labels = dom_lab) +
    geom_vline(data = block_bounds, aes(xintercept = xmax), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = mech_sep, linetype = "solid", color = "grey20", linewidth = 0.5) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25, linetype = "dashed"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.title = element_blank())
  p <- p + scale_size_area(max_size = 8, labels = scales::percent)
  if(square){
    p <- p + scale_fill_manual(values = mech_colors)
  } else {
    p <- p + scale_color_manual(values = mech_colors)
  }
  return(p)
}
p_dot <- plot_gene_domain_dot(genemeta); print(p_dot)
p_sq  <- plot_gene_domain_dot(genemeta, square = TRUE); print(p_sq)

######7 CADD and frequency
library(openxlsx); library(ggplot2); library(dplyr); library(patchwork)
genemeta <- read.xlsx("C:/Users/Administrator/Desktop/15个基因已发表突变中是JIA 表型的搜集/score_meta.xlsx")
data1<-genemeta[which(!is.na(genemeta$CADD)),]%>%filter(!is.na(GnomadV4))
plot_cadd_gnomad <- function(data, x_var = "CADD", point_size = 3, point_stroke = 0,
                             gene_order = c("TNFRSF1A","NLRP3","COPA","TNFAIP3","PSTPIP1","ADA2","LACC1","MVK","LPIN2","ELF4"),
                             title = "Summary of genes associated with arthritis"){
  library(ggplot2); library(dplyr)
  data$Gene <- factor(data$Gene, levels = gene_order)
  gene_colors <- c("TNFRSF1A"="#FEC260","NLRP3"="#FEC260","COPA"="#FEC260","TNFAIP3"="#FEC260","PSTPIP1"="#FEC260",
                   "ADA2"="#6699CC","LACC1"="#6699CC","MVK"="#6699CC","LPIN2"="#6699CC",
                   "ELF4"="#AF478A")
  if(x_var == "CADD"){
    p <- ggplot(data, aes(x = CADD, y = GnomadV4, fill = Gene)) +
      geom_point(size = point_size, shape = 21, stroke = point_stroke, color = "black") +
      geom_vline(xintercept = 20, linetype = "dashed", color = "grey50") +
      scale_y_log10() + scale_fill_manual(values = gene_colors, guide = "none") +
      labs(x = "CADD score", y = "gnomAD AF (log10)")
  } else {
    p <- ggplot(data, aes(x = GnomadV4, y = CADD, fill = Gene)) +
      geom_point(size = point_size, shape = 21, stroke = point_stroke, color = "black") +
      geom_hline(yintercept = 20, linetype = "dashed", color = "grey50") +
      scale_x_log10() + scale_fill_manual(values = gene_colors, guide = "none") +
      labs(x = "gnomAD AF", y = "CADD score")
  }
  p + facet_wrap(~ Gene, nrow = 2, scales = "free") +
    ggtitle(title) + theme_bw() +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid = element_blank(), legend.position = "none")
}
p1 <- plot_cadd_gnomad(data1, x_var = "CADD", point_size = 2.6,point_stroke = 0.5)
#p2 <- plot_cadd_gnomad(data1, x_var = "AF", point_stroke = 0.5)
p1
