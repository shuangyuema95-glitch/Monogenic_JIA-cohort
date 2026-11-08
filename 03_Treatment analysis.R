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
plot_drug_bar <- function(data, drug_class = NULL, bar_width = 0.7, color = "#6699CC") {
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

res <- plot_drug_bar(eff_tab, bar_width = 0.6, color = "#6699CC")
res$vertical      
res$horizontal   
res$pie           
res$facet_v       
res$facet_h

######2 sig use (no result was detected with p.adjust < 0.05)
library(openxlsx)
library(dplyr)

###(1) running analysis
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", check.names = FALSE)
colnames(treat) <- make.names(colnames(treat))
out_col <- grep("Response", colnames(treat), value = TRUE)[1]
eff_col <- grep("Effective", colnames(treat), value = TRUE)[1]
stages <- c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up")
drug_cols <- make.names(c("Antibiotics", "Non steroidal Anti-Inflammatory Drugs (NSAIDs)",
                          "Glucocorticoids", "Conventional synthetic DMARDs (csDMARDs)", "Anti-TNF", "Anti-IL-1",
                          "Anti-IL-6", "B-cell targeted therapy", "JAK inhibitors"))
drug_labs <- c("Antibiotics", "NSAIDs", "Glucocorticoids", "csDMARDs", "Anti-TNF", "Anti-IL-1", "Anti-IL-6", "B-cell", "JAKi")

run_fisher <- function(stage, outcome = c("response", "sfr"), gc_strict = FALSE) {
  outcome <- match.arg(outcome)
  data1 <- treat[treat$Stage == stage, ]
  stopifnot(nrow(data1) > 0)
  resp <- trimws(as.character(data1[[out_col]]))
  if (outcome == "response") {
    ev <- ifelse(resp %in% c("Complete response", "Partial response"), 1L,
                 ifelse(resp == "Nonresponse", 0L, NA_integer_))
  } else {
    gc <- trimws(as.character(data1[["Glucocorticoids"]]))
    gc_known <- !is.na(data1[["Glucocorticoids"]]) & gc != ""
    gu <- ifelse(gc_known, gc != "unused", NA)
    if (gc_strict) {
      e <- grepl("prednisolone|methylprednisolone|dexamethasone|prednisone|hydrocortisone",
                 tolower(as.character(data1[[eff_col]])))
      gu <- ifelse(gc_known, gu | e, NA)
    }
    resp_known <- resp %in% c("Complete response", "Partial response", "Nonresponse")
    ev <- ifelse(!gc_known | !resp_known, NA_integer_,
                 as.integer(!gu & resp %in% c("Complete response", "Partial response")))
  }
  ok <- !is.na(ev)
  res_list <- list()
  for (j in seq_along(drug_cols)) {
    if (outcome == "sfr" && drug_cols[j] == "Glucocorticoids") next
    x <- trimws(as.character(data1[[drug_cols[j]]]))
    used <- !is.na(x) & x != "unused"
    unused <- !is.na(x) & x == "unused"
    n_use <- sum(used & ok); n_unused <- sum(unused & ok)
    k_use <- sum(ev[used & ok] == 1); k_unused <- sum(ev[unused & ok] == 1)
    if (n_use == 0 || n_unused == 0) {
      res_list[[drug_labs[j]]] <- data.frame(Stage = stage, Outcome = outcome, Drug = drug_labs[j],
                                             n_total = sum(ok), used_n = n_use, event_used = k_use, non_event_used = n_use - k_use,
                                             event_used_pct = NA, unused_n = n_unused, event_unused = k_unused,
                                             non_event_unused = n_unused - k_unused, event_unused_pct = NA,
                                             OR = NA, p = NA, n_missing_drug = sum(ok & is.na(x)), n_missing_outcome = sum(!ok))
    } else {
      ft <- fisher.test(matrix(c(k_use, n_use - k_use, k_unused, n_unused - k_unused), nrow = 2, byrow = TRUE))
      res_list[[drug_labs[j]]] <- data.frame(Stage = stage, Outcome = outcome, Drug = drug_labs[j],
                                             n_total = sum(ok), used_n = n_use, event_used = k_use, non_event_used = n_use - k_use,
                                             event_used_pct = round(100 * k_use / n_use, 1), unused_n = n_unused, event_unused = k_unused,
                                             non_event_unused = n_unused - k_unused, event_unused_pct = round(100 * k_unused / n_unused, 1),
                                             OR = unname(ft$estimate), p = ft$p.value,
                                             n_missing_drug = sum(ok & is.na(x)), n_missing_outcome = sum(!ok))
    }
  }
  df <- do.call(rbind, res_list)
  df$p_adj <- p.adjust(df$p, method = "BH")
  df
}

res_resp_all <- bind_rows(lapply(stages, function(s) run_fisher(s, "response")))
res_sfr_all  <- bind_rows(lapply(stages, function(s) run_fisher(s, "sfr")))
res_resp_all
res_sfr_all
all(res_sfr_all$used_n + res_sfr_all$unused_n + res_sfr_all$n_missing_drug == res_sfr_all$n_total)
unique(res_sfr_all$n_total)
all(res_resp_all$used_n + res_resp_all$unused_n + res_resp_all$n_missing_drug == res_resp_all$n_total)
unique(res_resp_all$n_total)

screen_sig <- function(df, cut = 0.05, adj = TRUE) {
  col <- if (adj) "p_adj" else "p"
  out <- df[!is.na(df[[col]]) & df[[col]] < cut, ]
  out <- out[order(out[[col]]), ]
  out[, c("Stage", "Outcome", "Drug", "used_n", "event_used", "event_used_pct",
          "unused_n", "event_unused", "event_unused_pct", "OR", "p", "p_adj")]
}

sig_resp_adj <- screen_sig(res_resp_all, 0.05, TRUE);  sig_resp_adj
sig_resp_raw <- screen_sig(res_resp_all, 0.05, FALSE); sig_resp_raw
sig_sfr_adj  <- screen_sig(res_sfr_all, 0.05, TRUE);   sig_sfr_adj
sig_sfr_raw  <- screen_sig(res_sfr_all, 0.05, FALSE);  sig_sfr_raw

###(2) pheatmap or bar for visualization
plot_fisher_res <- function(df, stage, drug, fill_cols = c("#D8A0C0", "#6699CC"),
                            low_col = "#6699CC", high_col = "#DF5734", text_size = 3.2) {
  r <- df[df$Stage == stage & df$Drug == drug, ]
  if (nrow(r) != 1) stop("stage / drug not found or duplicated")
  ylab <- ifelse(r$Outcome == "sfr", "Steroid-free remission (%)", "Response (%)")
  plab <- paste0("Fisher p = ", formatC(r$p, format = "e", digits = 2))
  bar <- data.frame(Group = factor(c("Yes", "No"), levels = c("Yes", "No")),
                    n = c(r$used_n, r$unused_n), e = c(r$event_used, r$event_unused),
                    pct = c(r$event_used_pct, r$event_unused_pct))
  bar$lab <- paste0(bar$e, "/", bar$n)
  p_bar <- ggplot(bar, aes(Group, pct, fill = Group)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = lab), vjust = -0.4, size = text_size) +
    scale_fill_manual(values = fill_cols) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.14))) +
    labs(x = paste0(drug, " (Yes / No)"), y = ylab,
         title = paste0(stage, "\n", plab)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 9.5),
          axis.text = element_text(color = "black", size = 9),
          axis.title = element_text(color = "black", size = 9),
          legend.position = "none")
  tile <- data.frame(Group = factor(c("Yes", "Yes", "No", "No"), levels = c("No", "Yes")),
                     Out = factor(c("SFR", "no SFR", "SFR", "no SFR"), levels = c("SFR", "no SFR")),
                     n = c(r$event_used, r$used_n - r$event_used,
                           r$event_unused, r$unused_n - r$event_unused))
  p_tile <- ggplot(tile, aes(Out, Group, fill = n)) +
    geom_tile(color = "white", linewidth = 0.8, width = 0.95, height = 0.9) +
    geom_text(aes(label = n), size = text_size + 1) +
    scale_fill_gradient(low = low_col, high = high_col, name = "n") +
    labs(x = NULL, y = NULL, title = paste0(drug, " - ", stage, "\n", plab)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 9.5),
          axis.text = element_text(color = "black", size = 9),
          panel.grid = element_blank())
  list(bar = p_bar, tile = p_tile, row = r)
}

r_il1 <- plot_fisher_res(res_sfr_all, "Max response", "Anti-IL-1")
r_il1$bar
r_il1$tile
r_il1$row





#View(treat)
# filter(treat,Stage=="Maximum disease activity")->data1
# colnames(data1)[c(5:16,19)]
# index<-c(5:13)#drug types and treatment response
# ###check (both drugs and "unused")
# for(i in index){
#   print(table(data1[,i]))
# }

# run_fisher <- function(stage) {
#   data1 <- filter(treat, Stage == stage)
#   data1 <- data1[data1[[20]] %in% c("Complete response", "Partial response", "Nonresponse"), ]
#   index <- c(5:13); out_col <- 20
#   res_list <- list()
#   for(i in index) {
#     noall <- data1[which(data1[,i] == "unused"), out_col]
#     l1 <- sum(noall == "Complete response"); l2 <- sum(noall == "Partial response"); l3 <- sum(noall == "Nonresponse"); l0 <- l1 + l2
#     yesall <- data1[which(data1[,i] != "unused"), out_col]
#     k1 <- sum(yesall == "Complete response"); k2 <- sum(yesall == "Partial response"); k3 <- sum(yesall == "Nonresponse"); k0 <- k1 + k2
#     if (length(yesall) == 0 || length(noall) == 0) {
#       res_list[[colnames(data1)[i]]] <- data.frame(used_n = length(yesall), unused_n = length(noall), resp_used = k0, nonresp_used = k3, resp_unused = l0, nonresp_unused = l3, resp_used_pct = NA, resp_unused_pct = NA, OR = NA, p = NA)
#     } else {
#       tab <- matrix(c(k0, k3, l0, l3), nrow = 2, byrow = TRUE,
#                     dimnames = list(Use = c("used", "unused"), Outcome = c("Response", "Nonresponse")))
#       ft <- fisher.test(tab, alternative = "two.sided")
#       res_list[[colnames(data1)[i]]] <- data.frame(used_n = length(yesall), unused_n = length(noall), resp_used = k0, nonresp_used = k3, resp_unused = l0, nonresp_unused = l3, resp_used_pct = round(100 * k0 / length(yesall), 1), resp_unused_pct = round(100 * l0 / length(noall), 1), OR = ft$estimate, p = ft$p.value)
#     }
#   }
#   fisher_df <- do.call(rbind, res_list)
#   fisher_df$p_adj <- p.adjust(fisher_df$p, method = "BH")
#   fisher_df
# }
# unique(treat$Stage)
# res_max <- run_fisher("Maximum disease activity")  
# res_max
# res_maxres <- run_fisher("Max response") 
# res_maxres
# res_last <- run_fisher("Last follow-up") 
# res_last
# res_before <- run_fisher("Before genetic diagnosis") 
# res_before

######3 bubble plot for therapies and gene as well as pathway (efficacy rate= (CR+0.5*PR)/)
### JIA gene x drug class x efficacy (circle + square, normal + swapped)
library(openxlsx); library(dplyr); library(tidyr); library(ggplot2)

treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", check.names = FALSE)
metainfo <- read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")

library(openxlsx); library(dplyr); library(tidyr); library(ggplot2); library(grDevices)

treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", check.names = FALSE)
colnames(treat) <- make.names(colnames(treat))
metainfo <- read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")
colnames(metainfo) <- make.names(colnames(metainfo))

plot_bubble <- function(stage = "Maximum disease activity",
                        label = TRUE,
                        drug_cols = c("Anti-TNF", "Anti-IL-1", "Anti-IL-6", "JAK inhibitors",
                                      "B-cell targeted therapy",
                                      "Conventional synthetic DMARDs (csDMARDs)",
                                      "Glucocorticoids",
                                      "Non steroidal Anti-Inflammatory Drugs (NSAIDs)"),
                        drug_labs = c("Anti-TNF", "Anti-IL-1", "Anti-IL-6", "JAK inhibitors",
                                      "B-cell", "csDMARDs", "Glucocorticoids", "NSAIDs"),
                        group_order = c("Inflammasome", "NFKB", "Interferon regulation",
                                        "Immune metabolism", "IECD", "Uncategoried"),
                        fill_cols = NULL,
                        low_col = "#6699CC", high_col = "#DF5734",
                        size_range = c(3, 12), text_size = 2.5) {
  drug_cols <- make.names(drug_cols)
  out_col <- grep("Response", colnames(treat), value = TRUE)[1]
  data1 <- treat[treat$Stage == stage, ]
  data1$Number <- trimws(as.character(data1$Number))
  data1$gene <- toupper(trimws(as.character(metainfo$gene[match(data1$Number, metainfo$ID)])))
  data1$pathway <- as.character(metainfo$pathway[match(data1$Number, metainfo$ID)])
  long <- data1 %>%
    filter(!is.na(.data[[out_col]]) & .data[[out_col]] != "unknown") %>%
    pivot_longer(all_of(drug_cols), names_to = "class", values_to = "drug") %>%
    mutate(class = factor(drug_labs[match(class, drug_cols)], levels = drug_labs)) %>%
    filter(!is.na(drug) & trimws(drug) != "unused" & !is.na(gene) & !is.na(pathway))
  bubble_df <- long %>%
    group_by(gene, pathway, class) %>%
    summarise(total = n(),
              eff = (sum(.data[[out_col]] == "Complete response") +
                       0.5 * sum(.data[[out_col]] == "Partial response")) / n(),
              .groups = "drop") %>%
    mutate(eff_pct = round(eff * 100, 1)) %>%
    filter(total > 0)
  bubble_df$pathway <- factor(bubble_df$pathway, levels = intersect(group_order, unique(bubble_df$pathway)))
  bubble_df <- bubble_df %>%
    group_by(pathway) %>%
    mutate(gene = factor(gene, levels = rev(sort(unique(gene))))) %>%
    ungroup()
  
  th <- theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = "gray90"),
              panel.grid.major.y = element_line(linetype = "dashed"),
              panel.grid.major.x = element_line(linetype = "dotted"),
              panel.border = element_rect(color = "black", linewidth = 0.5),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", hjust = 0, size = 7),
              axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
              axis.text.y = element_text(color = "black", face = "italic"),
              axis.ticks = element_line(color = "black"))
  base <- list(facet_grid(pathway ~ ., scales = "free_y", space = "free_y"),
               labs(x = NULL, y = NULL), theme_bw(), th)
  
  fs <- if (is.null(fill_cols)) scale_fill_gradient(low = low_col, high = high_col, limits = c(0, 100))
  else scale_fill_gradientn(colours = fill_cols, limits = c(0, 100))
  size_n <- scale_size_continuous(range = size_range, name = "Total Treated (n)",
                                  breaks = sort(unique(bubble_df$total)))
  
  p_circle <- ggplot(bubble_df, aes(x = class, y = gene)) +
    geom_point(aes(size = total, fill = eff_pct), shape = 21, color = "black", stroke = 0.5) +
    fs + size_n + base
  if (label) p_circle <- p_circle + geom_text(aes(label = eff_pct), size = text_size)
  
  p_square <- ggplot(bubble_df, aes(x = class, y = gene)) +
    geom_point(aes(size = total, fill = eff_pct), shape = 22, color = "black", stroke = 0.5) +
    fs + size_n + base
  if (label) p_square <- p_square + geom_text(aes(label = eff_pct), size = text_size)
  
  tot_rng <- range(bubble_df$total)
  if (tot_rng[1] == tot_rng[2]) tot_rng <- tot_rng + c(-0.5, 0.5)
  fs_rev <- if (is.null(fill_cols)) scale_fill_gradient(low = low_col, high = high_col, limits = tot_rng)
  else scale_fill_gradientn(colours = fill_cols, limits = tot_rng)
  size_eff <- scale_size_continuous(range = size_range, name = "Efficacy Rate (%)",
                                    breaks = sort(unique(bubble_df$eff_pct)))
  
  p_circle_rev <- ggplot(bubble_df, aes(x = class, y = gene)) +
    geom_point(aes(size = eff_pct, fill = total), shape = 21, color = "black", stroke = 0.5) +
    fs_rev + size_eff + base
  if (label) p_circle_rev <- p_circle_rev + geom_text(aes(label = total), size = text_size)
  
  p_square_rev <- ggplot(bubble_df, aes(x = class, y = gene)) +
    geom_point(aes(size = eff_pct, fill = total), shape = 22, color = "black", stroke = 0.5) +
    fs_rev + size_eff + base
  if (label) p_square_rev <- p_square_rev + geom_text(aes(label = total), size = text_size)
  
  list(circle = p_circle, square = p_square,
       circle_rev = p_circle_rev, square_rev = p_square_rev,
       data = bubble_df)
}

cols_list <- list(p1 = colorRampPalette(c("#6699CC", "#DF5734"))(101))

stages <- c("Before genetic diagnosis", 
            "Maximum disease activity", 
            "Max response", "Last follow-up")
res_all <- list()
for (stg in stages) {
  res_all[[stg]] <- plot_bubble(stage = stg, fill_cols = cols_list$p1, label = TRUE)
}

(res_all[["Maximum disease activity"]]$circle+theme(legend.position = "none"))|res_all[["Max response"]]$circle


######4 upset plots for combined therapies and line plot for disease  activity score 
### JIA combination therapy module: UpSet + gene-labeled score + pathway heat/pie remove non-response  only including partial and complete response
library(openxlsx); library(dplyr); library(UpSetR); library(ggplot2)
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", check.names = FALSE)
colnames(treat) <- make.names(colnames(treat))
clin <- read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")

class_cols <- c(5:13)
short_names <- c("Antibiotics", "NSAIDs", "Glucocorticoids", "csDMARDs",
                 "Anti-TNF", "Anti-IL-1", "Anti-IL-6", "B-cell", "JAKi")
treat2 <- treat
colnames(treat2)[class_cols] <- short_names
score_cols <- 21
pathway_col <- c("NFKB" = "#B3D1E7", "Inflammasome" = "#6699CC",
                 "Immune metabolism" = "#847AB3", "Interferon regulation" = "#EC706E",
                 "IECD" = "#FEC260", "Uncategoried" = "#EBB1A4")

prep_combo <- function(stage, score_col = score_cols) {
  data1 <- treat2[treat2$Stage == stage, ]
  out_col <- grep("Response", colnames(data1), value = TRUE)[1]
  data1 <- data1[trimws(as.character(data1[[out_col]])) %in% c("Complete response", "Partial response"), ]
  bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
  colnames(bin) <- short_names
  keep <- rowSums(bin) >= 2
  bin <- bin[keep, , drop = FALSE]
  score <- as.numeric(data1[[score_col]][keep])
  id <- trimws(as.character(data1$Number[keep]))
  pathway <- clin$pathway[match(id, clin$ID)]
  gene <- toupper(trimws(as.character(clin$gene[match(id, clin$ID)])))
  ok <- !is.na(pathway) & !is.na(score)
  list(bin = bin[ok, , drop = FALSE], score = score[ok],
       pathway = pathway[ok], gene = gene[ok], id = id[ok])
}

combo_pw_long <- function(d, nintersects = 10) {
  set_order <- names(sort(colSums(d$bin), decreasing = TRUE))
  set_order <- set_order[colSums(d$bin[, set_order, drop = FALSE]) > 0]
  n <- length(set_order)
  main_list <- list(); meta_list <- list()
  for (code in 1:(2^n - 1)) {
    bits <- as.integer(intToBits(code))[1:n]
    sel <- set_order[rev(bits) == 1]
    idx <- rowSums(d$bin[, sel, drop = FALSE]) == length(sel) &
      rowSums(d$bin[, setdiff(short_names, sel), drop = FALSE]) == 0
    f <- sum(idx)
    if (f > 0) {
      pw_n <- table(d$pathway[idx])
      for (p in names(pw_n)) main_list[[length(main_list) + 1]] <- data.frame(
        combo = paste(sel, collapse = " & "), pathway = p, n = as.numeric(pw_n[p]))
      meta_list[[length(meta_list) + 1]] <- data.frame(
        combo = paste(sel, collapse = " & "), total = f,
        genes = paste(sort(unique(d$gene[idx])), collapse = ", "),
        patients = paste(sort(d$id[idx]), collapse = ", "),
        mean_score = mean(d$score[idx], na.rm = TRUE))
    }
  }
  meta_df <- do.call(rbind, meta_list)
  meta_df <- meta_df[order(meta_df$total, decreasing = TRUE), ]
  keep <- head(meta_df$combo, nintersects)
  long <- do.call(rbind, main_list)
  long <- long[long$combo %in% keep, ]
  long$combo <- factor(long$combo, levels = keep)
  long$pathway <- factor(long$pathway, levels = names(pathway_col))
  list(long = long, meta = meta_df[meta_df$combo %in% keep, ])
}

plot_combo_pathway_pie <- function(stage, nintersects = 10, pie_cols = pathway_col,
                                   ncol_pie = 5, show_pct = TRUE, pct_size = 2.2,
                                   label_r = 1.55, strip_size = 7,
                                   legend_size = 7, legend_key = 0.3,
                                   score_col = score_cols) {
  d <- prep_combo(stage, score_col)
  l <- combo_pw_long(d, nintersects)
  df <- as.data.frame(l$long)
  df <- df %>% group_by(combo) %>%
    arrange(desc(as.character(pathway))) %>%
    mutate(pct = n / sum(n) * 100, ypos = cumsum(pct) - pct / 2) %>%
    ungroup()
  p <- ggplot(df, aes(x = 1, y = pct, fill = pathway)) +
    geom_col(width = 1) +
    coord_polar(theta = "y", clip = "off") +
    facet_wrap(~combo, ncol = ncol_pie) +
    scale_fill_manual(values = pie_cols, name = NULL) +
    theme_void() +
    theme(strip.text = element_text(size = strip_size),
          legend.position = "bottom", legend.key.size = unit(legend_key, "cm"),
          legend.text = element_text(size = legend_size))
  if (show_pct) p <- p + geom_text(aes(x = label_r, y = ypos, label = paste0(round(pct, 1), "%")),
                                   size = pct_size, inherit.aes = FALSE)
  list(plot = p, data = df)
}

plot_combo_module <- function(stage, nintersects = 10,
                              main_color = "#6699CC", sets_color = "#847AB3",
                              line_color = "#AF478A", point_size = 2, gene_text_size = 2.2,
                              heat_col = "#6699CC", heat_na = "gray90", heat_text_size = 2.5,
                              pie_cols = pathway_col, ncol_pie = 5,
                              pie_show_pct = TRUE, pie_pct_size = 2.2,
                              pie_label_r = 1.55, pie_strip_size = 7,
                              pie_legend_size = 7, pie_legend_key = 0.3,
                              score_col = score_cols) {
  d <- prep_combo(stage, score_col)
  l <- combo_pw_long(d, nintersects)
  combo_df <- l$meta
  combo_df$combo <- factor(combo_df$combo, levels = as.character(combo_df$combo))
  
  ymax <- ceiling(max(combo_df$mean_score, na.rm = TRUE) / 2) * 2
  p_score <- ggplot(combo_df, aes(x = combo, y = mean_score, group = 1)) +
    geom_hline(yintercept = c(3, 7), linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(color = line_color, linewidth = 0.8) +
    geom_point(color = line_color, size = point_size) +
    geom_text(aes(label = genes), vjust = -1.5, size = gene_text_size) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, by = 2),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Mean score") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"))
  
  set_order <- names(sort(colSums(d$bin), decreasing = TRUE))
  set_order <- set_order[colSums(d$bin[, set_order, drop = FALSE]) > 0]
  p_upset <- upset(as.data.frame(d$bin), sets = set_order, nsets = length(set_order),
                   nintersects = nintersects, order.by = "freq",
                   main.bar.color = main_color, sets.bar.color = sets_color,
                   mainbar.y.label = "No. of patients", sets.x.label = "Patients using class")
  
  grid <- expand.grid(combo = levels(l$long$combo), pathway = names(pathway_col))
  dfm <- merge(grid, as.data.frame(l$long), by = c("combo", "pathway"), all.x = TRUE)
  dfm$combo <- factor(dfm$combo, levels = levels(l$long$combo))
  dfm$pathway <- factor(dfm$pathway, levels = rev(names(pathway_col)))
  dfm$fill_lab <- ifelse(is.na(dfm$n), "absent", "present")
  p_heat <- ggplot(dfm, aes(x = combo, y = pathway, fill = fill_lab)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(n), "", n)), size = heat_text_size) +
    scale_fill_manual(values = c(present = heat_col, absent = heat_na), guide = "none") +
    labs(x = NULL, y = "Pathway") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks = element_blank())
  
  p_pie <- plot_combo_pathway_pie(stage, nintersects, pie_cols = pie_cols, ncol_pie = ncol_pie,
                                  show_pct = pie_show_pct, pct_size = pie_pct_size,
                                  label_r = pie_label_r, strip_size = pie_strip_size,
                                  legend_size = pie_legend_size, legend_key = pie_legend_key,
                                  score_col = score_col)$plot
  
  table_out <- combo_df %>%
    transmute(Combo = as.character(combo), N = total, Mean_score = round(mean_score, 2),
              Genes = genes, Patients = patients)
  list(upset = p_upset, score = p_score, heat = p_heat, pie = p_pie, table = table_out)
}

stages <- c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up")
res_combo <- lapply(stages, function(s) plot_combo_module(
  stage = s, nintersects = 10,
  main_color = "#6699CC", sets_color = "#847AB3", line_color = "#AF478A", point_size = 2,
  gene_text_size = 2.2, heat_col = "#6699CC", heat_na = "gray90", heat_text_size = 2.5,
  pie_cols = pathway_col, ncol_pie = 5, pie_show_pct = TRUE, pie_pct_size = 2.2,
  pie_label_r = 1.55, pie_strip_size = 7, pie_legend_size = 7, pie_legend_key = 0.3,
  score_col = 21))
names(res_combo) <- stages

library(patchwork); library(ggplotify)
res_combo[["Maximum disease activity"]]$upset
res_combo[["Max response"]]$upset
p2 <- res_combo[["Maximum disease activity"]]$score | res_combo[["Max response"]]$score
p3 <- res_combo[["Maximum disease activity"]]$heat  | res_combo[["Max response"]]$heat
p4 <- res_combo[["Maximum disease activity"]]$pie   | res_combo[["Max response"]]$pie

p1; p2; p3; p4
res_combo[["Maximum disease activity"]]$table
res_combo[["Max response"]]$table

#######5 Sangke plots for targeted therapy
library(openxlsx); library(dplyr); library(ggplot2); library(ggalluvial)

treat <- read.xlsx("E:/Cohort PPT/JIA/code/treatment data.xlsx")
clin  <- read.xlsx("E:/Cohort PPT/JIA/0903_Clean case.xlsx")

class_cols <- 5:13
short_names <- c("Antibiotics","NSAIDs","Glucocorticoids","csDMARDs","Anti-TNF","Anti-IL-1","Anti-IL-6","B-cell","JAKi")
treat2 <- treat
colnames(treat2)[class_cols] <- short_names

pathway_col <- c(NFKB = "#B3D1E7", Inflammasome = "#6699CC", `Immune metabolism` = "#847AB3",
                 `Interferon regulation` = "#EC706E", IECD = "#FEC260", Uncategoried = "#EBB1A4")
drug_tar <- c("Anti-IL-1", "Anti-IL-6", "Anti-TNF", "JAKi")
stages <- c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up")

sankey_df <- function(stage, drug_cols = drug_tar, min_class = 2,
                      drop_nonresponse = TRUE, exclude_unused = TRUE) {
  data1 <- filter(treat2, Stage == stage)
  bin <- sapply(class_cols, function(i) as.integer(!is.na(data1[[i]]) & data1[[i]] != "unused"))
  colnames(bin) <- short_names
  keep <- rowSums(bin) >= min_class
  mid <- apply(data1[keep, drug_cols, drop = FALSE], 1, function(r) {
    v <- trimws(r[!is.na(r) & trimws(r) != "unused"])
    if (length(v) == 0) "unused" else paste(v, collapse = " & ")
  })
  oc <- grep("Response", colnames(data1), value = TRUE)
  df <- data.frame(pathway = clin$pathway[match(data1$Number[keep], clin$ID)],
                   gene = clin$gene[match(data1$Number[keep], clin$ID)],
                   combo = mid,
                   outcome = trimws(as.character(data1[[oc]][keep])),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$pathway) & df$outcome != "unknown", ]
  if (drop_nonresponse) df <- df[df$outcome != "Nonresponse", ]
  if (exclude_unused) df <- df[df$combo != "unused", ]
  df
}

plot_sankey3 <- function(df, by = c("pathway", "gene"), width = 0.44,
                         colors = pathway_col, label_size = 2.5) {
  by <- match.arg(by)
  df$A1 <- factor(df[[by]])
  df$combo <- factor(df$combo)
  df$outcome <- factor(df$outcome)
  p <- ggplot(df, aes(axis1 = A1, axis2 = combo, axis3 = outcome)) +
    geom_alluvium(aes(fill = pathway), width = width) +
    geom_stratum(width = width, fill = "grey92", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = label_size) +
    scale_x_discrete(limits = c(ifelse(by == "pathway", "Pathway", "Gene"), "Drugs", "Outcome"),
                     expand = c(0.05, 0.05)) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none")
  list(plot = p, data = df)
}

sk_df <- sankey_df("Max response")
sk_pw <- plot_sankey3(sk_df, "pathway"); sk_pw$plot
sk_gn <- plot_sankey3(sk_df, "gene"); sk_gn$plot

plot_sankey4 <- function(df, width = 0.36, colors = pathway_col, label_size = 2.0) {
  df$pathway <- factor(df$pathway)
  df$gene <- factor(df$gene)
  df$combo <- factor(df$combo)
  df$outcome <- factor(df$outcome)
  p <- ggplot(df, aes(axis1 = pathway, axis2 = gene, axis3 = combo, axis4 = outcome)) +
    geom_alluvium(aes(fill = pathway), width = width) +
    geom_stratum(width = width, fill = "grey92", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = label_size) +
    scale_x_discrete(limits = c("Pathway", "Gene", "Drugs", "Outcome"), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none")
  list(plot = p, data = df)
}

sk_4 <- plot_sankey4(sk_df); sk_4$plot

sankey_bundle <- function(stage, min_class = 2) {
  df <- sankey_df(stage, min_class = min_class)
  list(pathway = plot_sankey3(df, "pathway"),
       gene = plot_sankey3(df, "gene"),
       four = plot_sankey4(df),
       data = df)
}

sk_all <- lapply(stages, sankey_bundle)
names(sk_all) <- stages
data.frame(Stage = stages, n = sapply(sk_all, function(x) nrow(x$data)))

sk_all[["Maximum disease activity"]]$pathway$plot
sk_all[["Maximum disease activity"]]$gene$plot
sk_all[["Maximum disease activity"]]$four$plot
sk_all[["Max response"]]$pathway$plot
sk_all[["Max response"]]$gene$plot
sk_all[["Max response"]]$four$plot