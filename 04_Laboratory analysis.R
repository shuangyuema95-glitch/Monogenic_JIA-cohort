library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

####1 bubble plot showing the ratio of abnormal laboratoty testing index for 2 stages, grouped by genes
library(openxlsx); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales)

lab_path  <- "E:/Cohort PPT/JIA/code/Labo.xlsx"
gene_path <- "E:/Cohort PPT/JIA/0903_Clean case.xlsx"

up_ind    <- c("CRP","ESR","SAA","IL-6","TNF-α","IL-1β","Total bilirubin","Fibrinogen","LDH")
down_ind  <- c("C3","Serum creatinine","Hemoglobin","Platelet")
stage_vec <- c("Maximum disease activity","Max response")

bar_col_up   <- c("#E66145","#FCBB84")
bar_col_down <- c("#6699CC","#A0B8E0")
pathway_col<-c("NFKB"="#B3D1E7","Inflammasome"="#6699CC",
               "Immune metabolism"="#847AB3","Interferon regulation"="#EC706E",
               "IECD"="#FEC260","Uncategoried"="#EBB1A4")
bubble_fill_pal <- list(up = c("#FDECE3","#F29A57","#BE1F1B"), down = c("#E9EAEC","#B3D1E7","#6699CC"))

geom_stroke       <- 0.3
bar_width_main    <- 0.5
bar_width_stack   <- 0.5
text_size_label   <- 2.7
base_font_size    <- 9
legend_key_h      <- unit(0.45,"cm")
legend_key_w      <- unit(0.3,"cm")

ref_tbl <- tribble(~ind,~lo,~hi,
                   "WBC",4,10,"Neutrophil count",2.5,7.5,"Hemoglobin",120,160,"Platelet",100,300,"Reticulocyte count",24,84,
                   "CRP",0.068,8.2,"ESR",0,20,"SAA",-Inf,10,"Ferritin",15,300,"Triglyceride",1.7,2.3,
                   "LDH",109,245,"Fibrinogen",2,4,"AST",-Inf,40,"ALT",-Inf,50,"Total bilirubin",0.1,1.2,
                   "IL-6",0,7,"IL-18",100,300,"TNF-α",0,10,"IL-1β",0,5,"CXCL10",0,7,"CCL3",5,50,
                   "C4",0.16,0.48,"C3",0.9,1.8,"IgG",7,16,"IgA",0.7,4.0,"IgM",0.6,2.5,"IgE",21,200,
                   "NK cell activity",47.6,76.8,"CD4/CD8",0.9,2.0,"Urine PCR",0,30,"Serum creatinine",100,110)

read_labo <- function(path){
  raw <- read.xlsx(path, colNames = FALSE)
  nm <- trimws(sub("\\(.*","",as.character(unlist(raw[2,]))))
  nm[nm=="Serum amyloid A"] <- "SAA"
  nm[nm=="Lymphocyte"] <- "CD4/CD8"
  nm[nm=="Urine protein/creatinine ratio"] <- "Urine PCR"
  colnames(raw) <- nm
  raw <- raw[-(1:2),]
  colnames(raw)[1:3] <- c("Number","Name","Stage")
  raw %>% mutate(Number = trimws(as.character(Number)), Stage = trimws(as.character(Stage)))
}

read_gene <- function(path){
  g <- read.xlsx(path, sheet = 1, colNames = TRUE) %>%
    select(ID, gene, pathway) %>%
    mutate(across(everything(), ~trimws(as.character(.x)))) %>%
    filter(!is.na(ID), ID != "", !is.na(gene), gene != "") %>%
    distinct(ID, .keep_all = TRUE)
  g
}

gene_level_cluster <- function(stat_df, gene_tbl, pw_order){
  gene_tbl <- gene_tbl %>% distinct(gene, pathway)
  gene_ord <- c()
  path_boundary_y <- c()
  for(p in pw_order){
    genes_p <- gene_tbl %>% filter(pathway == p) %>% pull(gene) %>% unique()
    if(length(genes_p) == 0) next
    sub <- stat_df %>% filter(gene %in% genes_p)
    if(nrow(sub) == 0){
      ord_p <- genes_p
    }else if(length(genes_p)>=2){
      m_wide <- sub %>% pivot_wider(id_cols = gene, names_from = c(ind,Stage), values_from = pct, values_fill = 0)
      mat <- m_wide %>% select(-gene) %>% as.matrix()
      hc <- hclust(dist(mat))
      ord_p <- m_wide$gene[hc$order]
    }else{
      ord_p <- genes_p
    }
    gene_ord <- c(gene_ord, ord_p)
    path_boundary_y <- c(path_boundary_y, length(gene_ord))
  }
  list(gene_order = gene_ord, boundary = path_boundary_y[-length(path_boundary_y)])
}

calc_bubble <- function(lab, g, ref_tbl, up_ind, down_ind, stage_vec){
  all_ind <- c(up_ind, down_ind)
  long <- lab %>%
    select(Number, Stage, any_of(ref_tbl$ind)) %>%
    pivot_longer(any_of(ref_tbl$ind), names_to = "ind", values_to = "val") %>%
    mutate(val = suppressWarnings(as.numeric(trimws(as.character(val))))) %>%
    filter(Stage %in% stage_vec, !is.na(val), ind %in% all_ind) %>%
    left_join(ref_tbl, by = "ind") %>%
    mutate(dir = if_else(ind %in% up_ind, "up", "down"), abn_high = as.integer(val > hi), abn_low  = as.integer(val < lo)) %>%
    inner_join(g, by = c("Number" = "ID"))
  
  stat_gene <- long %>%
    group_by(gene, pathway, ind, dir, Stage) %>%
    summarise(n = n(), pct = 100 * if(first(dir)=="up") sum(abn_high)/n() else sum(abn_low)/n(), .groups = "drop") %>%
    filter(pct > 0) %>%
    mutate(ind = factor(ind, levels = all_ind), Stage = factor(Stage, levels = stage_vec))
  
  stat_global <- long %>%
    group_by(ind, Stage) %>%
    summarise(n_total = n(), pct_high = round(100*sum(abn_high)/n_total,2), pct_low  = round(100*sum(abn_low)/n_total,2), .groups = "drop") %>%
    mutate(ind = factor(ind, levels = all_ind), Stage = factor(Stage, levels = stage_vec))
  
  stack_raw <- long %>%
    mutate(is_abn = if_else(dir=="up", abn_high, abn_low)) %>%
    filter(is_abn == 1) %>%
    group_by(ind, Stage, pathway) %>%
    summarise(cnt = n(), .groups = "drop")
  
  stack_dat <- stack_raw %>%
    complete(ind, Stage, pathway = names(pathway_col), fill = list(cnt=0)) %>%
    group_by(ind, Stage) %>%
    mutate(total_abn = sum(cnt), pct_path = ifelse(total_abn == 0, 0, round(100 * cnt / total_abn, 2))) %>%
    ungroup() %>%
    mutate(ind = factor(ind, levels = all_ind), Stage = factor(Stage, levels = stage_vec))
  
  list(gene = stat_gene, global = stat_global, stack = stack_dat, long_raw = long)
}

theme_bub <- function(base_size){
  theme_classic(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 8),
          axis.text.y = element_text(colour = "black", size = 8),
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0.02,"lines"), plot.margin = margin(2,2,2,2),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}

mk_panel <- function(d, inds, cols, gene_lv, path_boundary, show_y, ttl, facet, stroke_val, base_size){
  d$gene <- factor(d$gene, levels = gene_lv)
  p <- ggplot(d, aes(x = if (facet) Stage else ind, y = gene)) +
    geom_hline(yintercept = path_boundary + 0.5, colour = "grey70", linetype = 2, linewidth =0.3)
  if (facet){
    link_data <- d %>%
      group_by(gene, pathway, ind) %>%
      filter(n_distinct(Stage) == 2) %>%
      pivot_wider(id_cols = c(gene, pathway, ind), names_from = Stage, values_from = c(n,pct), names_sep = "_") %>%
      filter(`n_Maximum disease activity` > 1 & `n_Max response` > 1) %>%
      mutate(trend = ifelse(`pct_Max response` < `pct_Maximum disease activity`, "decrease", "increase"), linetype = ifelse(trend=="decrease","solid","dashed")) %>%
      pivot_longer(cols = starts_with(c("pct_","n_")), names_to = c(".value","Stage"), names_sep = "_") %>%
      mutate(Stage = factor(Stage, levels = stage_vec))
    p <- p + geom_line(data = link_data, aes(x = Stage, y = gene, group = interaction(gene, ind), linetype = linetype), color = "black", linewidth = 0.3) +
      scale_linetype_manual(values = c("solid"="solid","dashed"="dashed"), guide = "none")
  }
  p <- p + geom_point(aes(size = n, fill = pct), shape = 21, colour = "grey25", stroke = stroke_val) +
    scale_fill_gradientn(colours = cols, limits = c(0,100), name = "Abnormal (%)") +
    scale_size(range = c(1.6, 7), name = "n", breaks = c(1,2,3,4), labels = ~paste0("n=(",.x,")")) +
    scale_y_discrete(limits = gene_lv, drop = FALSE) +
    theme_bub(base_size) + labs(title = ttl)
  if (facet){
    p <- p + facet_wrap(~ind, nrow = 1, scales = "free_x", drop = FALSE) +
      scale_x_discrete(drop = FALSE) + theme(strip.background = element_blank(), strip.text = element_text(size = 8))
  }else{
    p <- p + scale_x_discrete(limits = inds, drop = FALSE)
  }
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  p
}

plot_stage <- function(dat, stage_lbl, up_ind, down_ind, gene_lv, path_bd, pal, stroke_val, base_size){
  d  <- dat %>% filter(Stage == .env$stage_lbl) %>% mutate(gene = as.character(gene), ind = as.character(ind))
  pu <- mk_panel(filter(d, ind %in% .env$up_ind),   up_ind,   pal$up,   gene_lv, path_bd, TRUE,  "Elevated", FALSE, stroke_val, base_size)
  pd <- mk_panel(filter(d, ind %in% .env$down_ind), down_ind, pal$down, gene_lv, path_bd, FALSE, "Decreased", FALSE, stroke_val, base_size)
  pu + pd + plot_layout(widths = c(length(up_ind), length(down_ind)), guides = "collect")
}

plot_comb <- function(dat, up_ind, down_ind, gene_lv, path_bd, pal, stroke_val, base_size){
  d  <- dat %>% mutate(gene = as.character(gene), ind = as.character(ind))
  pu <- mk_panel(filter(d, ind %in% .env$up_ind)   %>% mutate(ind = factor(ind, levels = .env$up_ind)),
                 up_ind, pal$up, gene_lv, path_bd, TRUE, "Elevated", TRUE, stroke_val, base_size)
  pd <- mk_panel(filter(d, ind %in% .env$down_ind) %>% mutate(ind = factor(ind, levels = .env$down_ind)),
                 down_ind, pal$down, gene_lv, path_bd, FALSE, "Decreased", TRUE, stroke_val, base_size)
  pu + pd + plot_layout(widths = c(length(up_ind), length(down_ind)), guides = "collect")
}

plot_bar_combined <- function(global_df, up_ind, down_ind, col_up, col_down, bar_w, txt_size, base_size){
  d_up <- global_df %>% filter(ind %in% up_ind) %>% mutate(ind = factor(ind, levels = up_ind))
  p1 <- ggplot(d_up, aes(x = ind, y = pct_high, fill = Stage)) +
    geom_col(position = position_dodge(width=0.8), width = bar_w) +
    geom_text(aes(label=paste0(sprintf("%.2f",pct_high),"\n(n=",n_total,")")), position=position_dodge(0.8), vjust=-0.3, size=txt_size) +
    scale_fill_manual(values = col_up) +
    scale_y_continuous(limits = c(0,100), expand = c(0,2)) +
    theme_classic(base_size = base_size) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8), legend.position = "none") +
    labs(y = "Abnormal (%)", x = NULL)
  
  d_down <- global_df %>% filter(ind %in% down_ind) %>% mutate(ind = factor(ind, levels = down_ind))
  p2 <- ggplot(d_down, aes(x = ind, y = pct_low, fill = Stage)) +
    geom_col(position = position_dodge(width=0.8), width = bar_w) +
    geom_text(aes(label=paste0(sprintf("%.2f",pct_low),"\n(n=",n_total,")")), position=position_dodge(0.8), vjust=-0.3, size=txt_size) +
    scale_fill_manual(values = col_down) +
    scale_y_continuous(limits = c(0,100), expand = c(0,2)) +
    theme_classic(base_size = base_size) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), legend.position = "none") +
    labs(x = NULL)
  p1 + p2 + plot_layout(widths = c(length(up_ind), length(down_ind)))
}

plot_stack_pathway <- function(stack_df, up_ind, down_ind, pw_col, bar_w_stack, base_size){
  d_up <- stack_df %>% filter(ind %in% up_ind) %>% mutate(ind = factor(ind, levels = up_ind))
  s1 <- ggplot(d_up, aes(x = Stage, y = pct_path/100, fill = pathway)) +
    geom_col(position = position_fill(), width=bar_w_stack, na.rm=TRUE) +
    scale_fill_manual(values = pw_col) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)) +
    facet_wrap(~ind, nrow = 1, strip.position = "bottom") +
    theme_classic(base_size = base_size) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8), strip.placement = "outside", strip.background = element_blank(), legend.position = "bottom") +
    labs(y = "Pathway % among abnormal cases", x = NULL)
  
  d_down <- stack_df %>% filter(ind %in% down_ind) %>% mutate(ind = factor(ind, levels = down_ind))
  s2 <- ggplot(d_down, aes(x = Stage, y = pct_path/100, fill = pathway)) +
    geom_col(position = position_fill(), width=bar_w_stack, na.rm=TRUE) +
    scale_fill_manual(values = pw_col) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)) +
    facet_wrap(~ind, nrow = 1, strip.position = "bottom") +
    theme_classic(base_size = base_size) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8), strip.placement = "outside", strip.background = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), legend.position = "none") +
    labs(x = NULL)
  s1 + s2 + plot_layout(widths = c(length(up_ind), length(down_ind)))
}

lab  <- read_labo(lab_path)
gen  <- read_gene(gene_path)
res_list <- calc_bubble(lab, gen, ref_tbl, up_ind, down_ind, stage_vec)

stat_gene  <- res_list$gene
stat_global <- res_list$global
stat_stack  <- res_list$stack
long_raw    <- res_list$long_raw

clust_res <- gene_level_cluster(stat_gene, gen, names(pathway_col))
glv <- clust_res$gene_order
path_bd <- clust_res$boundary

p_act  <- plot_stage(stat_gene, stage_vec[1], up_ind, down_ind, glv, path_bd, bubble_fill_pal, geom_stroke, base_font_size)
p_res  <- plot_stage(stat_gene, stage_vec[2], up_ind, down_ind, glv, path_bd, bubble_fill_pal, geom_stroke, base_font_size)
p_comb <- plot_comb(stat_gene, up_ind, down_ind, glv, path_bd, bubble_fill_pal, geom_stroke, base_font_size)
p_bar_all <- plot_bar_combined(stat_global, up_ind, down_ind, bar_col_up, bar_col_down, bar_width_main, text_size_label, base_font_size)
p_stack_pathway <- plot_stack_pathway(stat_stack, up_ind, down_ind, pathway_col, bar_width_stack, base_font_size)

anno_theme <- theme(legend.position = "right", legend.key.height = legend_key_h, legend.key.width = legend_key_w,
                    legend.title = element_text(size=8), legend.text = element_text(size=7))

p_act <- p_act + plot_annotation(theme = anno_theme)
p_res <- p_res + plot_annotation(theme = anno_theme)
p_comb <- p_comb + plot_annotation(theme = anno_theme)

p_act
p_res
p_comb
p_bar_all
p_stack_pathway

stat_result<-list(
  stat_gene=stat_gene,
  stat_global=stat_global,
  stat_stack=stat_stack,
  long_raw=long_raw)

write.xlsx(stat_result,file="E:\\Cohort PPT\\JIA\\code\\1113_Labbubblebar_stats.xlsx")


######2 correlation analysis for clinical signatures and laboratory findings
library(openxlsx)
library(dplyr)
library(tidyr)
library(psych)
library(pheatmap)

cor_heat_merge_wrapper <- function(
    lab_path,
    clin_path,
    stage,
    lab_ref,
    drop_clin = c("Complement system dysregulation/activation",
                  "SLE overlap / Systemic lupus erythematosus",
                  "Positive autoantibodies","Crohn’s disease",
                  "Neutropenia","Macrophage activation syndrome (MAS)","Vitamin D deficiency"),
    sig_thr = 0.05,
    colors = c("#542788", "#F7F7F7", "#D6604D"),
    transpose = FALSE,
    border_col = "white",
    fontsize = 7
){
  read_labo <- function(path){
    labo <- read.xlsx(path, colNames = FALSE)
    ind <- as.character(unlist(labo[2, ]))
    colnames(labo) <- trimws(sub("\\(.*", "", ind))
    labo <- labo[-(1:2), ]
    colnames(labo)[1:3] <- c("Number", "Name", "Stage")
    colnames(labo)[colnames(labo) == "Serum amyloid A"] <- "SAA"
    colnames(labo)[colnames(labo) == "Lymphocyte"] <- "CD4/CD8"
    colnames(labo)[colnames(labo) == "Urine protein/creatinine ratio"] <- "Urine PCR"
    labo$Stage <- trimws(as.character(labo$Stage))
    labo
  }
  
  lab_pick <- function(labo, lab_ref,
                       stages = c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up")) {
    cls <- names(lab_ref)
    num_inds <- unlist(lapply(cls, function(cl) if (is.list(lab_ref[[cl]])) names(lab_ref[[cl]]) else NULL))
    bin_inds <- unlist(lapply(cls, function(cl) if (!is.list(lab_ref[[cl]])) lab_ref[[cl]] else NULL))
    pick <- function(x4) {
      ok <- !is.na(x4)
      if (sum(ok) == 0) return(NA)
      if (!is.na(x4["Maximum disease activity"])) return(x4["Maximum disease activity"])
      if (sum(ok) == 1) return(x4[ok])
      max(x4, na.rm = TRUE)
    }
    num_code <- function(x) {
      x <- trimws(as.character(x))
      x[x == "not available"] <- NA
      suppressWarnings(as.numeric(x))
    }
    bin_code <- function(x) { v <- trimws(tolower(as.character(x))); ifelse(v == "positive", 1, ifelse(v == "negative", 0, NA)) }
    mk <- function(inds, code) {
      out <- sapply(inds, function(cc) {
        w <- labo %>% select(Number, Stage, val = all_of(cc))
        w$val <- code(w$val)
        w <- w %>% pivot_wider(names_from = Stage, values_from = val)
        apply(w[, stages, drop = FALSE], 1, pick)
      })
      rownames(out) <- unique(labo$Number)
      out
    }
    list(numeric = mk(num_inds, num_code), binary = mk(bin_inds, bin_code))
  }
  
  lab_pick_binary <- function(labo, lab_ref, stage,
                              dec_inds = c("C3", "C4", "NK cell activity"),
                              both_inds = c("IgG", "IgM", "IgA")) {
    num_ref <- bind_rows(lapply(names(lab_ref), function(cl) {
      if (!is.list(lab_ref[[cl]])) return(NULL)
      bind_rows(lapply(names(lab_ref[[cl]]), function(m) {
        data.frame(ind = m, lo = lab_ref[[cl]][[m]][1], hi = lab_ref[[cl]][[m]][2],
                   stringsAsFactors = FALSE)
      }))
    }))
    lo <- setNames(num_ref$lo, num_ref$ind)
    hi <- setNames(num_ref$hi, num_ref$ind)
    num_code <- function(x) {
      x <- trimws(as.character(x))
      x[x == "not available"] <- NA
      suppressWarnings(as.numeric(x))
    }
    sub <- labo[labo$Stage == stage, ]
    all_inds <- c(dec_inds, both_inds, setdiff(num_ref$ind, c(dec_inds, both_inds)))
    cols <- unlist(lapply(all_inds, function(cc) {
      if (cc %in% both_inds) paste0(c("elevated", "decreased"), " ", cc)
      else if (cc %in% dec_inds) paste0("decreased ", cc)
      else paste0("elevated ", cc)
    }))
    m <- matrix(NA_real_, nrow = nrow(sub), ncol = length(cols),
                dimnames = list(sub$Number, cols))
    for (j in seq_along(all_inds)) {
      cc <- all_inds[j]
      if (!cc %in% colnames(sub)) next
      v <- num_code(sub[[cc]])
      if (cc %in% dec_inds) {
        m[, j] <- ifelse(!is.na(v), ifelse(v < lo[cc], 1, 0), NA)
      } else if (cc %in% both_inds) {
        m[, j] <- ifelse(!is.na(v), ifelse(v > hi[cc], 1, 0), NA)
        m[, match(paste0("decreased ", cc), colnames(m))] <-
          ifelse(!is.na(v), ifelse(v < lo[cc], 1, 0), NA)
      } else {
        m[, j] <- ifelse(!is.na(v), ifelse(v > hi[cc], 1, 0), NA)
      }
    }
    keep <- apply(m, 2, function(x)
      sum(x == 1, na.rm = TRUE) > 0 & sum(x == 0, na.rm = TRUE) > 0)
    m[, keep, drop = FALSE]
  }
  
  corr_binary <- function(matA, matB) {
    common <- intersect(rownames(matA), rownames(matB))
    matA <- matA[common, , drop = FALSE]
    matB <- matB[common, , drop = FALSE]
    out <- list(); k <- 0
    for (a in colnames(matA)) for (b in colnames(matB)) {
      x <- matA[, a]; y <- matB[, b]
      ok <- !is.na(x) & !is.na(y)
      n <- sum(ok)
      if (n < 3) next
      ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      if (!is.finite(ct$estimate)) next
      k <- k + 1
      out[[k]] <- data.frame(var1 = a, var2 = b, n = n,
                             rho = as.numeric(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE)
    }
    if (k == 0) return(data.frame(var1 = character(), var2 = character(), n = integer(),
                                  rho = numeric(), p = numeric(), stringsAsFactors = FALSE))
    res <- bind_rows(out)
    res$p_adj <- p.adjust(res$p, method = "BH")
    res[order(res$p_adj, res$p), ]
  }
  
  labo <- read_labo(lab_path)
  mat <- lab_pick(labo, lab_ref)
  mat_bin_ab <- mat$binary
  
  raw <- read.xlsx(clin_path, colNames = FALSE)
  pats <- as.character(raw[1, 3:ncol(raw)])
  feats <- as.character(raw[14:85, 2])
  m <- raw[14:85, 3:ncol(raw)]
  m <- as.matrix(m)
  code <- function(x) { v <- trimws(tolower(as.character(x))); ifelse(v == "yes", 1, ifelse(v == "no", 0, NA)) }
  clin_mat <- t(matrix(as.numeric(code(m)), nrow = nrow(m),
                       dimnames = list(feats, pats)))
  keep <- colSums(!is.na(clin_mat)) > 0
  clin_mat <- clin_mat[, keep, drop = FALSE]
  clin_mat <- clin_mat[, order(colSums(clin_mat, na.rm = TRUE), decreasing = TRUE), drop = FALSE]
  keep_clin <- !trimws(colnames(clin_mat)) %in% drop_clin
  clin_mat <- clin_mat[, keep_clin, drop = FALSE]
  
  stages <- c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up")
  mat_num_up <- setNames(lapply(stages, function(st) lab_pick_binary(labo, lab_ref, st)), stages)
  
  m1 <- mat_bin_ab
  m2 <- mat_num_up[[stage]]
  common_id <- intersect(intersect(rownames(m1), rownames(m2)), rownames(clin_mat))
  m1 <- m1[common_id, , drop = FALSE]
  m2 <- m2[common_id, , drop = FALSE]
  clin_sub <- clin_mat[common_id, , drop = FALSE]
  x_merge <- cbind(m1, m2)
  cor_res <- corr.test(x = x_merge, y = clin_sub, method = "spearman", use = "pairwise", adjust = "none")
  r <- cor_res$r; p <- cor_res$p
  sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  keep_r <- rowSums(p < sig_thr, na.rm = TRUE) > 0
  keep_c <- colSums(p < sig_thr, na.rm = TRUE) > 0
  if (any(keep_r) & any(keep_c)) { r <- r[keep_r, keep_c, drop = FALSE]; p <- p[keep_r, keep_c, drop = FALSE]; sig <- sig[keep_r, keep_c, drop = FALSE] }
  sig[is.na(p)] <- ""
  
  # filter: remove rows/columns with ONLY ONE significant result
  sig_bool <- p < sig_thr
  sig_bool[is.na(sig_bool)] <- FALSE
  row_sig_cnt <- rowSums(sig_bool)
  col_sig_cnt <- colSums(sig_bool)
  keep_r2 <- row_sig_cnt > 1
  keep_c2 <- col_sig_cnt > 1
  r <- r[keep_r2, keep_c2, drop = FALSE]
  p <- p[keep_r2, keep_c2, drop = FALSE]
  sig <- sig[keep_r2, keep_c2, drop = FALSE]
  
  max_abs <- max(abs(r), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  if (transpose) {
    r <- t(r)
    sig <- t(sig)
  }
  row_group <- ifelse(rownames(r) %in% colnames(m1), "Autoantibody", "Lab_abnormal")
  ann_row <- data.frame(group = row_group, row.names = rownames(r))
  ann_col <- list(group = c(Autoantibody="#999999", Lab_abnormal="#333333"))
  pheatmap(r, color = colorRampPalette(colors)(100),
           breaks = seq(-max_abs, max_abs, length.out = 100),
           cluster_rows = FALSE, cluster_cols = TRUE,
           treeheight_row = 0, treeheight_col = 10,
           display_numbers = sig, number_color = "black",
           fontsize = fontsize, border_color = border_col,
           annotation_row = ann_row, annotation_colors = ann_col,
           angle_col = 45)
  invisible(list(r = r, p = p, sig = sig))
}

lab_ref <- list(
  "Complete Blood Count (CBC)" = list(WBC = c(4, 10), "Neutrophil count" = c(2.5, 7.5), Hemoglobin = c(120, 160), Platelet = c(100, 300), "Reticulocyte count" = c(24, 84)),
  "Inflammatory markers" = list(CRP = c(0.068, 8.2), ESR = c(0, 20), SAA = c(-Inf, 10)),
  "Autoantibody" = c("RF", "anti-CCP", "HLA-B27", "ANA", "dsDNA", "ANCA", "Anticardiolipin antibody"),
  "Metabolic" = list(Ferritin = c(15, 300), Triglyceride = c(1.7, 2.3), LDH = c(109, 245), Fibrinogen = c(2, 4)),
  "Liver function" = list(AST = c(-Inf, 40), ALT = c(-Inf, 50), "Total bilirubin" = c(0.1, 1.2)),
  "Cytokines" = list("IL-6" = c(0, 7), "IL-18" = c(100, 300), "TNF-α" = c(0, 10), "IL-1β" = c(0, 5)),
  "Chemokines" = list(CXCL10 = c(0, 7), CCL3 = c(5, 50)),
  "Complement" = list(C4 = c(0.16, 0.48), C3 = c(0.9, 1.8)),
  "Immunoglobulin" = list(IgG = c(7, 16), IgA = c(0.7, 4.0), IgM = c(0.6, 2.5), IgE = c(21, 200)),
  "Immune function" = list("NK cell activity" = c(47.6, 76.8), "CD4/CD8" = c(0.9, 2.0)),
  "Renal function" = list("Urine PCR" = c(0, 30), "Serum creatinine" = c(100, 110)))

res <- cor_heat_merge_wrapper(
  lab_path = "E:\\Cohort PPT\\JIA\\code\\Labo.xlsx",
  clin_path = "E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx",
  stage = "Maximum disease activity",
  lab_ref = lab_ref,
  sig_thr = 0.05,
  colors = c("#542788", "#F7F7F7", "#D6604D"),
  transpose = FALSE,
  border_col = "white",
  fontsize = 7
)

######3 treatment and lab tesing 
#' Fisher test for association between medication and lab improvement (MDA -> MR)
#' @param lab_path path of Labo.xlsx
#' @param treat_path path of treatment data.xlsx
#' @param st_mda label for maximum disease activity stage
#' @param st_mr label for max response stage
#' @param lab_ref nested list containing reference range for continuous indicators
#' @return tibble of fisher test results, including OR, p, BH adjusted p
assoc_drug_lab_improve <- function(
    lab_path,
    treat_path,
    st_mda = "Maximum disease activity",
    st_mr = "Max response",
    lab_ref
) {
  library(openxlsx); library(dplyr); library(tidyr); library(stringr)
  
  ref_flat <- unlist(lab_ref, recursive = FALSE)
  is_rng   <- vapply(ref_flat, function(x) is.numeric(x) && length(x) == 2, logical(1))
  ref_df <- data.frame(
    ind = sub("^.*\\.", "", names(ref_flat)[is_rng]),
    lo  = vapply(ref_flat[is_rng], `[`, numeric(1), 1),
    hi  = vapply(ref_flat[is_rng], `[`, numeric(1), 2),
    stringsAsFactors = FALSE)
  
  labo <- read.xlsx(lab_path, colNames = FALSE)
  nm <- trimws(sub("\\(.*", "", as.character(unlist(labo[2, ]))))
  nm <- make.unique(nm, sep = "_")
  nm[1:3] <- c("Number", "Name", "Stage")
  colnames(labo) <- nm
  labo <- labo[-(1:2), ]
  colnames(labo)[colnames(labo) == "Serum amyloid A"]            <- "SAA"
  colnames(labo)[colnames(labo) == "Lymphocyte"]                 <- "CD4/CD8"
  colnames(labo)[colnames(labo) == "Urine protein/creatinine ratio"] <- "Urine PCR"
  labo$Number <- trimws(as.character(labo$Number))
  labo$Stage  <- trimws(as.character(labo$Stage))
  
  exist   <- intersect(ref_df$ind, colnames(labo))
  ref_sub <- ref_df[match(exist, ref_df$ind), ]
  
  d1 <- labo[labo$Stage == st_mda, c("Number", exist), drop = FALSE]
  d2 <- labo[labo$Stage == st_mr,  c("Number", exist), drop = FALSE]
  d1 <- d1[!duplicated(d1$Number), , drop = FALSE]
  d2 <- d2[!duplicated(d2$Number), , drop = FALSE]
  for (cc in exist) {
    d1[[cc]] <- suppressWarnings(as.numeric(trimws(as.character(d1[[cc]]))))
    d2[[cc]] <- suppressWarnings(as.numeric(trimws(as.character(d2[[cc]]))))
  }
  
  pats <- union(d1$Number, d2$Number)
  imp_mat <- matrix(NA_real_, length(pats), length(exist), dimnames = list(pats, exist))
  for (cc in exist) {
    a  <- d1[[cc]][match(pats, d1$Number)]
    b  <- d2[[cc]][match(pats, d2$Number)]
    lo <- ref_sub$lo[ref_sub$ind == cc]; hi <- ref_sub$hi[ref_sub$ind == cc]
    v  <- rep(NA_real_, length(pats))
    ok <- !is.na(a) & !is.na(b)
    v[ok] <- 0
    v[ok & a > hi] <- as.numeric(b[ok & a > hi] < a[ok & a > hi])
    v[ok & a < lo] <- as.numeric(b[ok & a < lo] > a[ok & a < lo])
    imp_mat[, cc] <- v
  }
  
  trt <- read.xlsx(treat_path)
  colnames(trt) <- gsub("\\.", " ", colnames(trt))
  eff_i <- grep("Effective Medications", colnames(trt))[1]
  eff <- data.frame(Number = trimws(as.character(trt[[1]])),
                    eff = trimws(as.character(trt[[eff_i]])), stringsAsFactors = FALSE)
  miss <- is.na(eff$eff) | eff$eff == "" | tolower(eff$eff) %in% c("na", "not available", "unknown", "n/a")
  has_rec <- tapply(!miss, eff$Number, any)
  
  eff_ok <- eff[!miss, ]
  eff_ok <- eff_ok %>% separate_rows(eff, sep = ",") %>% mutate(eff = trimws(eff)) %>%
    filter(eff != "", !tolower(eff) %in% c("no", "none"))
  eff_ok$key <- tolower(gsub("\\s+", " ", eff_ok$eff))
  drug_key <- eff_ok %>% distinct(key, .keep_all = TRUE) %>% transmute(key, drug = eff)
  
  drug_mat <- matrix(0L, length(pats), nrow(drug_key), dimnames = list(pats, drug_key$drug))
  hit <- eff_ok %>% distinct(Number, key)
  ri <- match(hit$Number, pats); ci <- match(hit$key, drug_key$key); ok_hit <- !is.na(ri)
  drug_mat[cbind(ri[ok_hit], ci[ok_hit])] <- 1L
  hr <- has_rec[pats]; hr[is.na(hr)] <- FALSE
  drug_mat[!hr, ] <- NA_integer_
  
  common <- intersect(rownames(imp_mat), rownames(drug_mat))
  res <- list(); k <- 0
  for (d in colnames(drug_mat)) for (m in colnames(imp_mat)) {
    x <- drug_mat[common, d]; y <- imp_mat[common, m]
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 4) next
    tb <- table(factor(x[ok], levels = c(0, 1)), factor(y[ok], levels = c(0, 1)))
    if (any(rowSums(tb) == 0) || any(colSums(tb) == 0)) next
    ft <- fisher.test(tb); k <- k + 1
    res[[k]] <- data.frame(drug = d, ind = m, n = sum(tb),
                           used_imp = tb[2, 2], used_noimp = tb[2, 1], unused_imp = tb[1, 2], unused_noimp = tb[1, 1],
                           n_used = sum(tb[2, ]), pct_imp_used = round(100 * tb[2, 2] / sum(tb[2, ]), 1),
                           pct_imp_unused = round(100 * tb[1, 2] / sum(tb[1, ]), 1),
                           OR = unname(ft$estimate), p = ft$p.value, stringsAsFactors = FALSE)
  }
  fisher_res <- bind_rows(res) %>% mutate(p_adj = p.adjust(p, method = "BH")) %>% arrange(p)
  return(fisher_res)
}

lab_path   <- "E:/Cohort PPT/JIA/code/Labo.xlsx"
treat_path <- "E:/Cohort PPT/JIA/code/treatment data.xlsx"
lab_ref <- list(
  "Complete Blood Count (CBC)" = list(WBC = c(4, 10), "Neutrophil count" = c(2.5, 7.5), Hemoglobin = c(120, 160), Platelet = c(100, 300), "Reticulocyte count" = c(24, 84)),
  "Inflammatory markers" = list(CRP = c(0.068, 8.2), ESR = c(0, 20), SAA = c(-Inf, 10)),
  "Autoantibody" = c("RF", "anti-CCP", "HLA-B27", "ANA", "dsDNA", "ANCA", "Anticardiolipin antibody"),
  "Metabolic" = list(Ferritin = c(15, 300), Triglyceride = c(1.7, 2.3), LDH = c(109, 245), Fibrinogen = c(2, 4)),
  "Liver function" = list(AST = c(-Inf, 40), ALT = c(-Inf, 50), "Total bilirubin" = c(0.1, 1.2)),
  "Cytokines" = list("IL-6" = c(0, 7), "IL-18" = c(100, 300), "TNF-α" = c(0, 10), "IL-1β" = c(0, 5)),
  "Chemokines" = list(CXCL10 = c(0, 7), CCL3 = c(5, 50)),
  "Complement" = list(C4 = c(0.16, 0.48), C3 = c(0.9, 1.8)),
  "Immunoglobulin" = list(IgG = c(7, 16), IgA = c(0.7, 4.0), IgM = c(0.6, 2.5), IgE = c(21, 200)),
  "Immune function" = list("NK cell activity" = c(47.6, 76.8), "CD4/CD8" = c(0.9, 2.0)),
  "Renal function" = list("Urine PCR" = c(0, 30), "Serum creatinine" = c(100, 110)))

fisher_res <- assoc_drug_lab_improve(lab_path, treat_path, lab_ref = lab_ref)
fisher_res %>% filter(p < 0.05)




# ######3 Specific laboratory findings for 4 stages
# plot_lab_panel <- function(labo, treat, pts, vars, line_col="#2E75B6", line_size=1, point_col="#C00000", point_size=2, stage_levels=c("Before genetic diagnosis","Maximum disease activity","Max response","Last follow-up"), stage_short=c("Before genetic diagnosis","Disease flare","Maximal response","Last follow-up")){
#   labo <- as.data.frame(labo, stringsAsFactors = FALSE)
#   if (!("Number" %in% names(labo)) && identical(trimws(as.character(labo[2,1])), "Number")) { names(labo) <- as.character(labo[2, ]); labo <- labo[-(1:2), ]; rownames(labo) <- NULL }
#   for (v in c("Number","Stage")) labo[[v]] <- as.character(labo[[v]])
#   norm <- function(x) gsub("[^[:alnum:]]","",tolower(x))
#   nm <- names(labo)
#   alias <- c(SAA="Serum amyloid A (SAA, <10 mg/L)","CD4/CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)","CD4.CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)",UPCR="Urine protein/creatinine ratio (0-30 mg/g)","Urine PCR"="Urine protein/creatinine ratio (0-30 mg/g)")
#   idx <- match(vars, nm)
#   if (any(is.na(idx))) { i2 <- match(norm(vars[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i2 }
#   v2 <- ifelse(vars %in% names(alias), unname(alias[vars]), vars)
#   if (any(is.na(idx))) { i3 <- match(v2[is.na(idx)], nm); idx[is.na(idx)] <- i3
#   if (any(is.na(idx))) { i4 <- match(norm(v2[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i4 } }
#   miss <- vars[is.na(idx)]
#   if (length(miss)) stop("columns not found in labo: ", paste(miss, collapse=", "))
#   use <- nm[idx]
#   for (v in use) labo[[v]] <- suppressWarnings(as.numeric(as.character(labo[[v]])))
#   labo <- labo[labo$Number %in% pts & labo$Stage %in% stage_levels, c("Number","Stage", use)]
#   long <- pivot_longer(labo, all_of(use), names_to="Indicator", values_to="Value")
#   long$Stage <- factor(long$Stage, levels=stage_levels)
#   long$Indicator <- factor(long$Indicator, levels=use)
#   dm <- data.frame(full=c("Prednisolone","Methylprednisolone","Methotrexate (MTX)","Mycophenolate mofetil (MMF)","Sulfasalazine","Cyclophosphamide","Leflunomide","Hydroxychloroquine (HCQ)","Cyclosporine","Tacrolimus","Adalimumab","Etanercept","Infliximab","Yisaipu","Anakinra","Canakinumab","Firsekibart","Tocilizumab","Belimumab","Tofacitinib","Baricitinib","Thalidomide","IVIG","Colchicine","Growth hormone"), abbr=c("PSL","mPSL","MTX","MMF","SSZ","CYC","LEF","HCQ","CsA","TAC","ADA","ETA","IFX","YSP","ANA","CAN","FIR","TCZ","BLM","TOF","BARI","THD","IVIG","COL","GH"), cat=c("GC","GC","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","Anti-TNF","Anti-TNF","Anti-TNF","Anti-TNF","Anti-IL-1","Anti-IL-1","Anti-IL-1","Anti-IL-6","B-cell","JAK","JAK","Immuno","IVIG","Colch","GH"), stringsAsFactors=FALSE)
#   dm$key <- gsub("\\s+","",tolower(dm$full))
#   cat_order <- c("GC","csDMARD","Anti-TNF","Anti-IL-1","Anti-IL-6","B-cell","JAK","Immuno","IVIG","Colch","GH")
#   tt <- as.data.frame(treat, stringsAsFactors = FALSE)
#   for (v in c("Number","Stage")) tt[[v]] <- as.character(tt[[v]])
#   cn <- names(tt)
#   pick <- function(pat) cn[grepl(pat, cn)][1]
#   c_glu <- pick("^Glucocort"); c_dm <- pick("DMARD"); c_tnf <- pick("Anti.TNF"); c_il1 <- pick("Anti.IL.1"); c_il6 <- pick("Anti.IL.6"); c_b <- pick("B.cell"); c_jak <- pick("^JAK"); c_ivig <- pick("IVIG"); c_col <- pick("^Colch"); c_gh <- pick("Growth"); c_imm <- pick("^Immunomod"); c_oth <- pick("^Other.supplemental")
#   parse_treat <- function(pid, stg) {
#     r <- tt[tt$Number == pid & tt$Stage == stg, , drop=FALSE]
#     if (!nrow(r)) return(NULL)
#     toks <- character(0)
#     take <- function(col, yesval=NULL) {
#       if (is.na(col)) return(invisible(NULL))
#       x <- trimws(as.character(r[[col]][1]))
#       if (is.na(x) || x %in% c("", "unused", "unknown")) return(invisible(NULL))
#       if (!is.null(yesval)) { if (x == "yes") toks <<- c(toks, yesval); return(invisible(NULL)) }
#       if (x == "yes") return(invisible(NULL))
#       for (pt in trimws(strsplit(x, ",")[[1]])) {
#         k <- gsub("\\s+", "", tolower(pt))
#         i <- match(k, dm$key)
#         if (is.na(i)) i <- match(sub("\\(.*\\)$", "", k), dm$key)
#         if (!is.na(i)) toks <<- c(toks, dm$abbr[i])
#       }
#       invisible(NULL)
#     }
#     take(c_glu); take(c_dm); take(c_tnf); take(c_il1); take(c_il6); take(c_b); take(c_jak); take(c_ivig, "IVIG"); take(c_col, "Colchicine"); take(c_gh, "Growth hormone"); take(c_imm); take(c_oth)
#     if (!length(toks)) return(NULL)
#     ab <- unique(toks)
#     ord <- order(match(dm$cat[match(ab, dm$abbr)], cat_order))
#     ab <- ab[ord]
#     c(paste(ab, collapse="+"), paste(dm$full[match(ab, dm$abbr)], collapse="+"))
#   }
#   drug_rows <- lapply(pts, function(p) {
#     m <- sapply(stage_levels, function(s) { r <- parse_treat(p, s); if (is.null(r)) c(NA_character_, NA_character_) else r })
#     data.frame(Number=p, Stage=stage_levels, Drug_abbr=m[1, ], Drug_full=m[2, ], stringsAsFactors=FALSE)
#   })
#   drug_df <- bind_rows(drug_rows)
#   drug_wide <- pivot_wider(drug_df, id_cols=Number, names_from=Stage, values_from=Drug_abbr)
#   ind_lbl <- setNames(sub(" \\(", "\n(", use), use)
#   plots <- list()
#   for (p in pts) {
#     dat <- long[long$Number == p, ]
#     dd <- drug_df[drug_df$Number == p, ]
#     xl <- sapply(seq_along(stage_levels), function(i) { d <- dd$Drug_abbr[i]; if (!is.na(d)) paste0(stage_short[i], "\n", d) else stage_short[i] })
#     plots[[p]] <- ggplot(dat, aes(Stage, Value)) + geom_line(aes(group=1), color=line_col, linewidth=line_size) + geom_point(color=point_col, size=point_size) + facet_wrap(~ Indicator, ncol=1, scales="free_y", strip.position="left", labeller=labeller(Indicator=ind_lbl)) + scale_x_discrete(labels=xl) + labs(title=p) + theme_bw() + theme(axis.title=element_blank(), axis.text.x=element_text(size=5.5, lineheight=0.85), strip.text.y.left=element_text(size=7, angle=0, lineheight=0.9), strip.background=element_rect(fill="grey92", linewidth=0.2), plot.title=element_text(size=10, face="bold", hjust=0.5, margin=margin(b=2)), panel.spacing=unit(0.15, "lines"), panel.grid.minor=element_blank())
#   }
#   list(plots=plots, drugs=drug_wide, drug_long=drug_df)
# }
# labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", sheetIndex = 1, header = FALSE)
# treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", sheetIndex = 1)
# pts <- c("P1","P3","P5","P7","P8","P9","P15","P16","P18","P20","P21","P25","P26");length(pts)
# vars <- c("Serum creatinine (100-110 μmol/L)","AST (<40 U/L)","ALT (<50 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)","Serum amyloid A (SAA, <10 mg/L)")
# res <- plot_lab_panel(labo, treat, pts=pts, vars=vars, line_col="#2E75B6", line_size=1.2, point_col="#C00000", point_size=2.5)
# names(res$plots)
# res$drugs
# res$plots$P16
# length(res$drugs)
# 
# library(ggpubr)
# library(gridExtra)
# n_p <- length(res$plots)
# setwd("E:\\Cohort PPT\\JIA\\code")
# pdf("lab_panel_all.pdf", width = 13.5, height = 9.5, onefile = TRUE)
# for (i in seq_len(ceiling(n_p / 4))) {
#   idx <- ((i - 1) * 4 + 1):min(i * 4, n_p)
#   grid.arrange(grobs = res$plots[idx], ncol = 4, respect = TRUE)
# }
# dev.off()
# 
# 
# suppressMessages({library(xlsx);library(dplyr);library(tidyr);library(ggplot2)})
# plot_lab_panel <- function(labo, treat, pts, vars, line_col="#2E75B6", line_size=1, point_size=2, col_high="red", col_low="blue", col_norm="black", stage_levels=c("Before genetic diagnosis","Maximum disease activity","Max response","Last follow-up"), stage_short=c("Before genetic diagnosis","Disease flare","Maximal response","Last follow-up")){
#   labo <- as.data.frame(labo, stringsAsFactors = FALSE)
#   if (!("Number" %in% names(labo)) && identical(trimws(as.character(labo[2,1])), "Number")) { names(labo) <- as.character(labo[2, ]); labo <- labo[-(1:2), ]; rownames(labo) <- NULL }
#   for (v in c("Number","Stage")) labo[[v]] <- as.character(labo[[v]])
#   norm <- function(x) gsub("[^[:alnum:]]","",tolower(x))
#   nm <- names(labo)
#   alias <- c(SAA="Serum amyloid A (SAA, <10 mg/L)","CD4/CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)","CD4.CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)",UPCR="Urine protein/creatinine ratio (0-30 mg/g)","Urine PCR"="Urine protein/creatinine ratio (0-30 mg/g)")
#   idx <- match(vars, nm)
#   if (any(is.na(idx))) { i2 <- match(norm(vars[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i2 }
#   v2 <- ifelse(vars %in% names(alias), unname(alias[vars]), vars)
#   if (any(is.na(idx))) { i3 <- match(v2[is.na(idx)], nm); idx[is.na(idx)] <- i3
#   if (any(is.na(idx))) { i4 <- match(norm(v2[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i4 } }
#   miss <- vars[is.na(idx)]
#   if (length(miss)) stop("columns not found in labo: ", paste(miss, collapse=", "))
#   use <- nm[idx]
#   for (v in use) labo[[v]] <- suppressWarnings(as.numeric(as.character(labo[[v]])))
#   labo <- labo[labo$Number %in% pts & labo$Stage %in% stage_levels, c("Number","Stage", use)]
#   long <- pivot_longer(labo, all_of(use), names_to="Indicator", values_to="Value")
#   long$Stage <- factor(long$Stage, levels=stage_levels)
#   long$Indicator <- factor(long$Indicator, levels=use)
#   ref <- data.frame(Indicator=c("Serum creatinine (100-110 μmol/L)","AST (<40 U/L)","ALT (<50 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)","Serum amyloid A (SAA, <10 mg/L)"), lo=c(100,0,0,0.068,0,0), hi=c(110,40,50,8.2,20,10), stringsAsFactors=FALSE)
#   long <- left_join(long, ref, by="Indicator")
#   long$lo[is.na(long$lo)] <- -Inf; long$hi[is.na(long$hi)] <- Inf
#   long$Status <- ifelse(long$Value > long$hi, "high", ifelse(long$Value < long$lo, "low", "normal"))
#   dm <- data.frame(full=c("Prednisolone","Methylprednisolone","Methotrexate (MTX)","Mycophenolate mofetil (MMF)","Sulfasalazine","Cyclophosphamide","Leflunomide","Hydroxychloroquine (HCQ)","Cyclosporine","Tacrolimus","Adalimumab","Etanercept","Infliximab","Yisaipu","Anakinra","Canakinumab","Firsekibart","Tocilizumab","Belimumab","Tofacitinib","Baricitinib","Thalidomide","IVIG","Colchicine","Growth hormone"), abbr=c("PSL","mPSL","MTX","MMF","SSZ","CYC","LEF","HCQ","CsA","TAC","ADA","ETA","IFX","YSP","ANA","CAN","FIR","TCZ","BLM","TOF","BARI","THD","IVIG","COL","GH"), cat=c("GC","GC","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","Anti-TNF","Anti-TNF","Anti-TNF","Anti-TNF","Anti-IL-1","Anti-IL-1","Anti-IL-1","Anti-IL-6","B-cell","JAK","JAK","Immuno","IVIG","Colch","GH"), stringsAsFactors=FALSE)
#   dm$key <- gsub("\\s+","",tolower(dm$full))
#   cat_order <- c("GC","csDMARD","Anti-TNF","Anti-IL-1","Anti-IL-6","B-cell","JAK","Immuno","IVIG","Colch","GH")
#   tt <- as.data.frame(treat, stringsAsFactors = FALSE)
#   for (v in c("Number","Stage")) tt[[v]] <- as.character(tt[[v]])
#   cn <- names(tt)
#   pick <- function(pat) cn[grepl(pat, cn)][1]
#   c_glu <- pick("^Glucocort"); c_dm <- pick("DMARD"); c_tnf <- pick("Anti.TNF"); c_il1 <- pick("Anti.IL.1"); c_il6 <- pick("Anti.IL.6"); c_b <- pick("B.cell"); c_jak <- pick("^JAK"); c_ivig <- pick("IVIG"); c_col <- pick("^Colch"); c_gh <- pick("Growth"); c_imm <- pick("^Immunomod"); c_oth <- pick("^Other.supplemental")
#   parse_treat <- function(pid, stg) {
#     r <- tt[tt$Number == pid & tt$Stage == stg, , drop=FALSE]
#     if (!nrow(r)) return(NULL)
#     toks <- character(0)
#     take <- function(col, yesval=NULL) {
#       if (is.na(col)) return(invisible(NULL))
#       x <- trimws(as.character(r[[col]][1]))
#       if (is.na(x) || x %in% c("", "unused", "unknown")) return(invisible(NULL))
#       if (!is.null(yesval)) { if (x == "yes") toks <<- c(toks, yesval); return(invisible(NULL)) }
#       if (x == "yes") return(invisible(NULL))
#       for (pt in trimws(strsplit(x, ",")[[1]])) {
#         k <- gsub("\\s+", "", tolower(pt))
#         i <- match(k, dm$key)
#         if (is.na(i)) i <- match(sub("\\(.*\\)$", "", k), dm$key)
#         if (!is.na(i)) toks <<- c(toks, dm$abbr[i])
#       }
#       invisible(NULL)
#     }
#     take(c_glu); take(c_dm); take(c_tnf); take(c_il1); take(c_il6); take(c_b); take(c_jak); take(c_ivig, "IVIG"); take(c_col, "Colchicine"); take(c_gh, "Growth hormone"); take(c_imm); take(c_oth)
#     if (!length(toks)) return(NULL)
#     ab <- unique(toks)
#     ord <- order(match(dm$cat[match(ab, dm$abbr)], cat_order))
#     ab <- ab[ord]
#     c(paste(ab, collapse="+"), paste(dm$full[match(ab, dm$abbr)], collapse="+"))
#   }
#   drug_rows <- lapply(pts, function(p) {
#     m <- sapply(stage_levels, function(s) { r <- parse_treat(p, s); if (is.null(r)) c(NA_character_, NA_character_) else r })
#     data.frame(Number=p, Stage=stage_levels, Drug_abbr=m[1, ], Drug_full=m[2, ], stringsAsFactors=FALSE)
#   })
#   drug_df <- bind_rows(drug_rows)
#   drug_wide <- pivot_wider(drug_df, id_cols=Number, names_from=Stage, values_from=Drug_abbr)
#   ind_lbl <- setNames(sub(" \\(", "\n(", use), use)
#   plots <- list()
#   for (p in pts) {
#     dat <- long[long$Number == p, ]
#     dd <- drug_df[drug_df$Number == p, ]
#     xl <- sapply(seq_along(stage_levels), function(i) { d <- dd$Drug_abbr[i]; if (!is.na(d)) paste0(stage_short[i], "\n", d) else stage_short[i] })
#     plots[[p]] <- ggplot(dat, aes(Stage, Value)) + geom_line(aes(group=1), color=line_col, linewidth=line_size) + geom_point(aes(color=Status), size=point_size) + facet_wrap(~ Indicator, ncol=1, scales="free_y", strip.position="top", labeller=labeller(Indicator=ind_lbl)) + scale_color_manual(values=c(high=col_high, normal=col_norm, low=col_low), name=NULL, breaks=c("high","normal","low"), labels=c("Above normal","Normal","Below normal")) + scale_x_discrete(labels=xl) + labs(title=p) + theme_classic() + theme(axis.title=element_blank(), axis.text.x=element_text(size=5.5, lineheight=0.85), axis.text.y=element_text(size=6), axis.line=element_line(color="grey40", linewidth=0.3), panel.grid.major.x=element_line(linetype="22", color="grey85", linewidth=0.4), panel.grid.major.y=element_line(linetype="11", color="grey85", linewidth=0.4), panel.grid.minor=element_blank(), strip.background=element_blank(), strip.text=element_text(face="bold", size=9), plot.title=element_text(size=11, face="bold", hjust=0.5, margin=margin(b=2)), panel.spacing=unit(0.2, "lines"), legend.position="top", legend.title=element_blank(), legend.text=element_text(size=7))
#   }
#   list(plots=plots, drugs=drug_wide, drug_long=drug_df)
# }
# labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", sheetIndex = 1, header = FALSE)
# treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", sheetIndex = 1)
# pts <- c("P5","P21")
# vars <- c("AST (<40 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)")
# 
# names(res$plots)
# 
# res <- plot_lab_panel(labo, treat, pts = c("P5","P21","P25"),
#                       vars = c("AST (<40 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)"),
#                       line_col = "gray60", line_size = 1.2,
#                       col_high = "#BB0021", 
#                       col_norm = "black", 
#                       col_low = "#6699CC",
#                       point_size = 3.5)
# library(patchwork)
# res$plots$P5|res$plots$P21


