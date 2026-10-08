###treatment analysis
## 8-10 non-response, 4-7 partial response; 0-3 complete response


######1 effective drug count (proportion )
library(ggplot2)
library(openxlsx)
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx")
eff <- treat$Effective.Medications.During.Treatment
all_effdrugs <- c()
for(i in 1:length(eff)){
  all_effdrugs <- c(all_effdrugs, trimws(unlist(strsplit(eff[i], split = ","))))
}
all_effdrugs <- all_effdrugs[!all_effdrugs %in% "no"]
sort(round(table(all_effdrugs)/(sum(table(all_effdrugs))),3),decreasing = T)


drug_class <- list("Antibiotics" = c("Penicillin", "Cephalosporins", "Carbapenems", "Macrolides", "Quinolones", "Glycopeptides"),
  "Non steroidal Anti-Inflammatory Drugs (NSAIDs)" = c("Ibuprofen", "Naproxen", "Diclofenac", "Indomethacin", "Celecoxib"),
  "Glucocorticoids" = c("Prednisolone", "Methylprednisolone", "Dexamethasone"),
  "Conventional synthetic DMARDs (csDMARDs)" = c("Methotrexate (MTX)", "Sulfasalazine", "Cyclosporine", "Azathioprine", "Cyclophosphamide", "Hydroxychloroquine (HCQ)", "Mycophenolate mofetil (MMF)"),
  "Anti-TNF" = c("Etanercept", "Adalimumab", "Infliximab", "Golimumab"),
  "Anti-IL-1" = c("Anakinra", "Canakinumab", "Firsekibart"),
  "Anti-IL-6" = c("Tocilizumab"),
  "B-cell targeted therapy" = c("Rituximab", "Belimumab"),
  "JAK inhibitors" = c("Tofacitinib", "Baricitinib", "Ruxolitinib", "Upadacitinib"),
  "Colchicine" = "Colchicine",
  "Immunomodulators" = "Thalidomide",
  "Others" = c("Intravenous immunoglobulin", "Growth hormone", "Vitamin D", "Furosemide", "Plasma exchange", "Antibiotics"))


map_df <- data.frame(drug = unlist(drug_class, use.names = FALSE), class = rep(names(drug_class), lengths(drug_class)))
eff_tab <- as.data.frame(table(all_effdrugs), stringsAsFactors = FALSE)
colnames(eff_tab) <- c("drug", "n")
eff_tab$class <- map_df$class[match(eff_tab$drug, map_df$drug)]
eff_tab$drug_prop <- round(eff_tab$n / sum(eff_tab$n), 4)
class_n <- tapply(eff_tab$n, eff_tab$class, sum)
eff_tab$class_prop <- round(class_n[eff_tab$class] / sum(eff_tab$n), 4)
eff_tab <- eff_tab[order(eff_tab$n, decreasing = TRUE), ]
eff_tab

library(ggplot2)
plot_drug_bar <- function(data, drug_class = NULL, bar_width = 0.7, color = "#847AB3") {
  if (is.data.frame(data)) df <- data.frame(drug = as.character(data[[1]]), n = as.numeric(data[[2]]), stringsAsFactors = FALSE)
  else if (!is.null(names(data))) df <- data.frame(drug = names(data), n = as.numeric(data), stringsAsFactors = FALSE)
  else { df <- as.data.frame(table(data), stringsAsFactors = FALSE); colnames(df) <- c("drug", "n") }
  df <- df[df$drug != "no" & df$n > 0, ]
  if (is.data.frame(data) && ncol(data) >= 3) df$class <- as.character(data[[3]])
  else if (!is.null(drug_class)) df$class <- rep(names(drug_class), lengths(drug_class))[match(df$drug, unlist(drug_class, use.names = FALSE))]
  if (any(is.na(df$class))) stop("unmapped drugs found")
  df$pct <- round(df$n / sum(df$n) * 100, 2)
  class_n <- tapply(df$n, df$class, sum)
  class_pct <- round(class_n / sum(df$n) * 100, 2)
  df$class_pct <- class_pct[df$class]
  
  df_d <- df[order(df$n, decreasing = TRUE), ]
  df_d$drug_v <- factor(df_d$drug, levels = df_d$drug)
  df_d$drug_h <- factor(df_d$drug, levels = rev(df_d$drug))
  p_v <- ggplot(df_d, aes(x = drug_v, y = pct)) +
    geom_col(width = bar_width, fill = color) +
    geom_text(aes(label = pct), vjust = -0.5, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"))
  p_h <- ggplot(df_d, aes(x = pct, y = drug_h)) +
    geom_col(width = bar_width, fill = color) +
    geom_text(aes(label = pct), hjust = -0.15, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Proportion (%)", y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"))
  
  class_df <- data.frame(class = names(class_n), n = as.numeric(class_n), stringsAsFactors = FALSE)
  class_df$pct <- round(class_df$n / sum(class_df$n) * 100, 2)
  class_df$class <- factor(class_df$class, levels = class_df$class[order(class_df$n, decreasing = TRUE)])
  p_pie <- ggplot(class_df, aes(x = "", y = n, fill = class)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = paste0(pct, "%")), position = position_stack(vjust = 0.5), size = 3) +
    coord_polar(theta = "y") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          legend.title = element_blank())
  
  df_f <- df[order(df$class_pct, df$n, decreasing = TRUE), ]
  df_f$class_f <- factor(df_f$class, levels = unique(df_f$class))
  df_f$drug_fv <- factor(df_f$drug, levels = unique(df_f$drug))
  df_f$drug_fh <- factor(df_f$drug, levels = unlist(lapply(unique(df_f$class), function(cl) rev(df_f$drug[df_f$class == cl])), use.names = FALSE))
  p_fv <- ggplot(df_f, aes(x = drug_fv, y = pct)) +
    geom_col(width = bar_width, fill = color) +
    geom_text(aes(label = pct), vjust = -0.5, size = 2.5) +
    facet_grid(. ~ class_f, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", size = 8))
  p_fh <- ggplot(df_f, aes(x = pct, y = drug_fh)) +
    geom_col(width = bar_width, fill = color) +
    geom_text(aes(label = pct), hjust = -0.15, size = 2.5) +
    facet_grid(class_f ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Proportion (%)", y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", size = 8))
  list(vertical = p_v, horizontal = p_h, pie = p_pie, facet_v = p_fv, facet_h = p_fh)
}

res <- plot_drug_bar(eff_tab, bar_width = 0.6, color = "#847AB3")
res$vertical      
res$horizontal   
res$pie           
res$facet_v       
res$facet_h

######2 sig use (no result was detected with p.adjust < 0.05)
library(openxlsx)
library(dplyr)
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx")

#View(treat)
# filter(treat,Stage=="Maximum disease activity")->data1
# colnames(data1)[c(5:16,19)]
# index<-c(5:13)#drug types and treatment response
# ###check (both drugs and "unused")
# for(i in index){
#   print(table(data1[,i]))
# }



run_fisher <- function(stage) {
  data1 <- filter(treat, Stage == stage)
  data1 <- data1[data1[[19]] %in% c("Complete response", "Partial response", "Nonresponse"), ]
  index <- c(5:13); out_col <- 19
  res_list <- list()
  for(i in index) {
    noall <- data1[which(data1[,i] == "unused"), out_col]
    l1 <- sum(noall == "Complete response"); l2 <- sum(noall == "Partial response"); l3 <- sum(noall == "Nonresponse"); l0 <- l1 + l2
    yesall <- data1[which(data1[,i] != "unused"), out_col]
    k1 <- sum(yesall == "Complete response"); k2 <- sum(yesall == "Partial response"); k3 <- sum(yesall == "Nonresponse"); k0 <- k1 + k2
    if (length(yesall) == 0 || length(noall) == 0) {
      res_list[[colnames(data1)[i]]] <- data.frame(used_n = length(yesall), unused_n = length(noall), resp_used = k0, nonresp_used = k3, resp_unused = l0, nonresp_unused = l3, resp_used_pct = NA, resp_unused_pct = NA, OR = NA, p = NA)
    } else {
      tab <- matrix(c(k0, k3, l0, l3), nrow = 2, byrow = TRUE,
                    dimnames = list(Use = c("used", "unused"), Outcome = c("Response", "Nonresponse")))
      ft <- fisher.test(tab, alternative = "two.sided")
      res_list[[colnames(data1)[i]]] <- data.frame(used_n = length(yesall), unused_n = length(noall), resp_used = k0, nonresp_used = k3, resp_unused = l0, nonresp_unused = l3, resp_used_pct = round(100 * k0 / length(yesall), 1), resp_unused_pct = round(100 * l0 / length(noall), 1), OR = ft$estimate, p = ft$p.value)
    }
  }
  fisher_df <- do.call(rbind, res_list)
  fisher_df$p_adj <- p.adjust(fisher_df$p, method = "BH")
  fisher_df
}
unique(treat$Stage)
res_max <- run_fisher("Maximum disease activity")  
res_max
res_maxres <- run_fisher("Max response") 
res_maxres
res_last <- run_fisher("Last follow-up") 
res_last
res_before <- run_fisher("Before genetic diagnosis") 
res_before


######3 upset plots for combined therapies and line plot for disease  activity score 
library(openxlsx)
library(dplyr)
library(UpSetR)
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx")
class_cols <- c(5:13)
colnames(treat)[class_cols]
short_names <- 
  c("Antibiotics", "NSAIDs", "Glucocorticoids", "csDMARDs",
    "Anti-TNF", "Anti-IL-1", "Anti-IL-6", "B-cell",
    "JAKi")
treat2 <- treat
colnames(treat2)[class_cols] <- short_names

plot_upset_stage <- function(stage, nintersects = 20,
                             main_color = "#6699CC", sets_color = "#847AB3") {
  data1 <- filter(treat2, Stage == stage)
  bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
  colnames(bin) <- short_names
  bin <- bin[rowSums(bin) >= 2, , drop = FALSE]
  cs <- sort(colSums(bin), decreasing = TRUE)
  set_order <- names(cs[cs > 0])
  upset(as.data.frame(bin), sets = set_order, nsets = length(set_order), nintersects = nintersects,
        order.by = "freq",
        main.bar.color = main_color, sets.bar.color = sets_color,
        mainbar.y.label = "No. of patients", sets.x.label = "Patients using class")
}

u1<-plot_upset_stage("Maximum disease activity")
u2<-plot_upset_stage("Max response")

score_cols<-grep("Disease.Activity.Score",colnames(treat2),value=t)
score_line <- function(stage, nintersects = 20, line_color = "#AF478A", score_col = score_cols) {
  data1 <- filter(treat2, Stage == stage)
  bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
  colnames(bin) <- short_names
  keep <- rowSums(bin) >= 2
  bin <- bin[keep, , drop = FALSE]
  score <- as.numeric(data1[[score_col]][keep])
  combo_list <- list()
  for (k in 2:length(short_names)) {
    com <- combn(short_names, k)
    for (j in seq_len(ncol(com))) {
      sel <- com[, j]
      idx <- rowSums(bin[, sel, drop = FALSE]) == k & rowSums(bin[, setdiff(short_names, sel), drop = FALSE]) == 0
      n <- sum(idx)
      if (n > 0) combo_list[[length(combo_list) + 1]] <- data.frame(combo = paste(sel, collapse = " & "), n = n, mean_score = mean(score[idx], na.rm = TRUE))
    }
  }
  combo_df <- do.call(rbind, combo_list)
  combo_df <- combo_df[order(combo_df$n, decreasing = TRUE), ]
  combo_df <- head(combo_df, nintersects)
  combo_df$combo <- factor(combo_df$combo, levels = combo_df$combo)
  ymax <- ceiling(max(combo_df$mean_score, na.rm = TRUE) / 2) * 2
  ggplot(combo_df, aes(x = combo, y = mean_score, group = 1)) +
    geom_hline(yintercept = c(3, 7), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(color = line_color, linewidth = 0.8) +
    geom_point(color = line_color, size = 2) +
    geom_text(aes(label = round(mean_score, 2)), vjust = -1, size = 3) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, by = 2), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Mean score") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"))
}
s1 <- score_line("Maximum disease activity")
s2 <- score_line("Max response")

###only one therapy (no retain)
# plot_upset_single <- function(stage, nintersects = 20, main_color = "#6699CC", sets_color = "#847AB3") {
#   data1 <- filter(treat2, Stage == stage)
#   bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
#   colnames(bin) <- short_names
#   bin <- bin[rowSums(bin) == 1, , drop = FALSE]
#   cs <- sort(colSums(bin), decreasing = TRUE)
#   set_order <- names(cs[cs > 0])
#   upset(as.data.frame(bin), sets = set_order, nsets = length(set_order), nintersects = nintersects,
#         order.by = "freq",
#         main.bar.color = main_color, sets.bar.color = sets_color,
#         mainbar.y.label = "No. of patients", sets.x.label = "Patients using class")
# }
# 
# score_line_single <- function(stage, nintersects = 20, line_color = "#AF478A", score_col = score_cols) {
#   data1 <- filter(treat2, Stage == stage)
#   bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
#   colnames(bin) <- short_names
#   keep <- rowSums(bin) == 1
#   bin <- bin[keep, , drop = FALSE]
#   score <- as.numeric(data1[[score_col]][keep])
#   combo_list <- list()
#   for (j in seq_along(short_names)) {
#     sel <- short_names[j]
#     idx <- bin[, sel] == 1
#     n <- sum(idx)
#     if (n > 0) combo_list[[length(combo_list) + 1]] <- data.frame(combo = sel, n = n, mean_score = mean(score[idx], na.rm = TRUE))
#   }
#   combo_df <- do.call(rbind, combo_list)
#   combo_df <- combo_df[order(combo_df$n, decreasing = TRUE), ]
#   combo_df <- head(combo_df, nintersects)
#   combo_df$combo <- factor(combo_df$combo, levels = combo_df$combo)
#   ymax <- ceiling(max(combo_df$mean_score, na.rm = TRUE) / 2) * 2
#   ggplot(combo_df, aes(x = combo, y = mean_score, group = 1)) +
#     geom_hline(yintercept = c(3, 7), linetype = "dashed", color = "grey50", linewidth = 0.5) +
#     geom_line(color = line_color, linewidth = 0.8) +
#     geom_point(color = line_color, size = 2) +
#     geom_text(aes(label = round(mean_score, 2)), vjust = -1, size = 3) +
#     scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, by = 2), expand = expansion(mult = c(0, 0.05))) +
#     labs(x = NULL, y = "Mean score") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
#           axis.text.y = element_text(color = "black", size = 7),
#           axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"))
# }
# 
# u1s <- plot_upset_single("Maximum disease activity")
# u2s <- plot_upset_single("Max response")
# s1s <- score_line_single("Maximum disease activity")
# s2s <- score_line_single("Max response")


#######4 Sangke plots for targeted therapy
library(dplyr)
library(ggplot2)
library(ggalluvial)

clin<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")

pathway_col <- c(
  "NFKB" = "#B3D1E7",
  "Inflammasome" = "#6699CC",
  "Immune metabolism" = "#847AB3",
  "Interferon" = "#EC706E",
  "Cell death" = "#FEC260",
  "Uncategoried" = "#EBB1A4")

sankey_drugcol <- function(stage, drug_cols, exclude_unused = FALSE,
                           width = 0.44, colors = pathway_col) {
  data1 <- filter(treat2, Stage == stage)
  bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
  colnames(bin) <- short_names
  keep <- rowSums(bin) >= 2
  mid <- apply(data1[keep, drug_cols, drop = FALSE], 1, function(r) {
    v <- trimws(r[!is.na(r) & trimws(r) != "unused"])
    if (length(v) == 0) "unused" else paste(v, collapse = " & ")
  })
  pathway <- clin$pathway[match(data1$Number[keep], clin$ID)]
  out_col <- grep("Response", colnames(data1), value = TRUE)
  outcome <- data1[[out_col]][keep]
  df <- data.frame(pathway = pathway, combo = mid, outcome = outcome, stringsAsFactors = FALSE)
  df <- df[!is.na(df$pathway) & df$outcome != "unknown", ]
  if (exclude_unused) df <- df[df$combo != "unused", ]
  df$pathway <- factor(df$pathway)
  df$combo <- factor(df$combo)
  df$outcome <- factor(df$outcome)
  lab2 <- df %>% count(combo, pathway) %>% group_by(combo) %>% mutate(pct = round(n / sum(n) * 100, 1), y = cumsum(n) - n / 2)
  lab3 <- df %>% count(outcome, pathway) %>% group_by(outcome) %>% mutate(pct = round(n / sum(n) * 100, 1), y = cumsum(n) - n / 2)
  lab4 <- df %>% count(outcome, combo) %>% group_by(outcome) %>% mutate(pct = round(n / sum(n) * 100, 1), y = cumsum(n) - n / 2)
  p <- ggplot(df, aes(axis1 = pathway, axis2 = combo, axis3 = outcome)) +
    geom_alluvium(aes(fill = pathway), width = width) +
    geom_stratum(width = width, fill = "grey92", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    geom_text(data = lab2, aes(x = 2, y = y, label = paste0(pct, "%")), inherit.aes = FALSE, size = 2.2) +
    geom_text(data = lab3, aes(x = 3, y = y, label = paste0(pct, "%")), inherit.aes = FALSE, size = 2.2) +
    scale_x_discrete(limits = c("Pathway", "Drugs", "Outcome"), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none")
  print(p)
  list(plot = p, drug_pathway_pct = lab2, outcome_pathway_pct = lab3, outcome_combo_pct = lab4)
}

res <- sankey_drugcol("Maximum disease activity",
                      c("Anti-IL-1", "Anti-IL-6", "Anti-TNF", "JAKi"), exclude_unused = TRUE)
res$plot
res$drug_pathway_pct
res$outcome_pathway_pct
res$outcome_combo_pct

