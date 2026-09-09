library(ggplot2)
library(openxlsx)
library(tidyverse)
metainfo<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")
country_df <- metainfo %>%
  count(country, name = "n") %>%
  mutate(percent = n / sum(n) * 100)

#####1 Figure 1  country proportion pie
country_colors <- c(
  "China"  = "gray90",
  "Egypt"  = "#B5BBE3",
  "Brazil" = "#CBDAA9")
p_country <- ggplot(country_df, aes(x = "", y = n, fill = country)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
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
sex_ration<-c(1,1,2,1,1,2,1,1,1,2,2,2,1,2,1,2,2,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,1)
sex_df<-data.frame(sex=sex_ration)|>count(sex,name="n")|>mutate(group=case_when(sex==2~"Female",sex==1~"Male"),percent=n/sum(n)*100)
sex_col<-c("Male"="#6699CC","Female"="#86C7B4")
p_sex<-ggplot(sex_df,aes(x="",y=n,fill=group))+geom_bar(stat="identity",width=1,color="white",linewidth=0.5)+coord_polar("y",start=0)+geom_text(aes(label=paste0(round(percent,1),"%")),position=position_stack(vjust=0.5),colour="black",size=4.2)+scale_fill_manual(values=sex_col)+theme_void()+theme(legend.title=element_blank(),legend.text=element_text(size=11,colour="black"),plot.title=element_text(hjust=0.5,size=14,colour="black"))
print(p_sex)

library(patchwork)
p_country|p_sex

######3 Bar plot for mutation type, novel or reported and genotype
library(tidyverse)
library(patchwork)

plot_meta_bar <- function(metainfo, target_col, fill_color = "#847AB3", bar_width = 0.7){
  split_vec <- metainfo[[target_col]] %>%
    str_split("/") %>%
    unlist() %>%
    str_trim()
  cnt <- tibble(val = split_vec) %>%
    count(val) %>%
    mutate(total = sum(n), percent = n / total * 100) %>%
    arrange(desc(percent))
  cnt$val <- factor(cnt$val, levels = cnt$val)
  
  p <- ggplot(cnt, aes(x = val, y = percent)) +
    geom_col(fill = fill_color, width = bar_width) +
    geom_text(aes(label = round(percent, 1)), vjust = -0.25, size = 3.5, colour = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
      axis.ticks.x = element_line(colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks.y = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(colour = "black")
    ) +
    labs(y = "Percentage (%)")
  return(p)
}

p_mut <- plot_meta_bar(metainfo, "MUT")
p_geno <- plot_meta_bar(metainfo, "genotype")+theme(axis.title.y = element_blank())
p_novel <- plot_meta_bar(metainfo, "novel")+theme(axis.title.y = element_blank())
p_acmg<-plot_meta_bar(metainfo,"ACMG_class")+theme(axis.title.y = element_blank())

p_mut|p_geno|p_novel|p_acmg


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


-####(2) density plots
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

