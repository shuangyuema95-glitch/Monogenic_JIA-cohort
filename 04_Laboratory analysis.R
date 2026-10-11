library(openxlsx)
library(dplyr)

######1 proportion for each laboratory testing
labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", colNames = FALSE)
labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", sheetIndex = 1, header = FALSE)
ind <- as.character(unlist(labo[2, ]))
colnames(labo) <- trimws(sub("\\(.*", "", ind))
labo <- labo[-(1:2), ]
colnames(labo)[1:3] <- c("Number", "Name", "Stage")
colnames(labo)[colnames(labo) == "Serum amyloid A"] <- "SAA"
colnames(labo)[colnames(labo) == "Lymphocyte"] <- "CD4/CD8"
colnames(labo)[colnames(labo) == "Urine protein/creatinine ratio"] <- "Urine PCR"

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

labo_prop <- function(labo, lab_ref, stages = NULL, digits = 1) {
  cls <- names(lab_ref)
  ref_num <- bind_rows(lapply(cls, function(cl) {
    if (!is.list(lab_ref[[cl]])) return(NULL)
    bind_rows(lapply(names(lab_ref[[cl]]), function(m) {
      lo <- lab_ref[[cl]][[m]][1]; hi <- lab_ref[[cl]][[m]][2]
      rng <- ifelse(is.infinite(lo) & lo < 0, paste0("<", hi),
                    ifelse(is.infinite(hi) & hi > 0, paste0(">", lo),
                           paste0(lo, "-", hi)))
      data.frame(class = cl, ind = m, lo = lo, hi = hi, ref_range = rng, stringsAsFactors = FALSE)
    }))
  }))
  ref_bin <- bind_rows(lapply(cls, function(cl) {
    if (is.list(lab_ref[[cl]])) return(NULL)
    data.frame(class = cl, ind = lab_ref[[cl]], stringsAsFactors = FALSE)
  }))
  drop_na <- function(x) trimws(as.character(x))[trimws(as.character(x)) != "not available"]
  calc <- function(sub) {
    a <- lapply(seq_len(nrow(ref_num)), function(i) {
      v <- suppressWarnings(as.numeric(drop_na(sub[[ref_num$ind[i]]])))
      n <- sum(!is.na(v))
      if (n == 0) return(NULL)
      data.frame(ind = ref_num$ind[i], class = ref_num$class[i], type = "numeric", n = n,
                 low = sum(v < ref_num$lo[i], na.rm = TRUE), high = sum(v > ref_num$hi[i], na.rm = TRUE),
                 normal = sum(v >= ref_num$lo[i] & v <= ref_num$hi[i], na.rm = TRUE),
                 pos = NA, neg = NA, ref_range = ref_num$ref_range[i], stringsAsFactors = FALSE)
    })
    b <- lapply(seq_len(nrow(ref_bin)), function(i) {
      v <- drop_na(sub[[ref_bin$ind[i]]])
      n <- sum(v %in% c("negative", "positive"))
      if (n == 0) return(NULL)
      data.frame(ind = ref_bin$ind[i], class = ref_bin$class[i], type = "binary", n = n,
                 low = NA, high = NA, normal = NA,
                 pos = sum(v == "positive"), neg = sum(v == "negative"),
                 ref_range = NA, stringsAsFactors = FALSE)
    })
    out <- c(a, b); out <- out[!vapply(out, is.null, logical(1))]
    bind_rows(out) %>% mutate(abn_pct = ifelse(type == "numeric", round(100 * (low + high) / n, digits), NA),
                              pos_pct = ifelse(type == "binary", round(100 * pos / n, digits), NA),
                              neg_pct = ifelse(type == "binary", round(100 * neg / n, digits), NA))
  }
  if (is.null(stages)) calc(labo) %>% mutate(stage = "all") %>% select(stage, everything())
  else bind_rows(lapply(stages, function(s) calc(filter(labo, Stage == s)) %>% mutate(stage = s))) %>% select(stage, everything())
}

res <- labo_prop(labo, lab_ref, stages = c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up"))

####check 
# colnames(labo)[4:11]
# res1<-data.frame()
# for(i in 4:11){
#   ll<-colnames(labo)[i]
#   range<-as.numeric(unlist(strsplit(res[match(ll,res$ind),"ref_range"],"-")))
#   quan<-labo[labo$Stage=="Before genetic diagnosis",ll]
#   quan<-quan[-grep("not available",quan)]
#   all<-length(quan)
#   quan<-as.numeric(quan)  
#   low<-sum(quan<range[1])
#   high<-sum(quan>range[2])
#   norm<-all-low-high
#   res1<-rbind(res1,c(ll,all,low,high,norm))
# }

library(ggplot2)
library(dplyr)

plot_lab_stack <- function(res, stage, colors = c(
  low = "#847AB3", high = "#EC706E", normal = "#B5BBE3",
  positive = "#AF478A", negative = "#6699CC"),
  bar_width = 0.87) {
  d <- res %>% filter(stage == !!stage)
  long <- bind_rows(
    d %>% filter(type == "numeric") %>% transmute(ind, class, n, segment = "low", pct = round(100 * low / n, 1)),
    d %>% filter(type == "numeric") %>% transmute(ind, class, n, segment = "normal", pct = round(100 * normal / n, 1)),
    d %>% filter(type == "numeric") %>% transmute(ind, class, n, segment = "high", pct = round(100 * high / n, 1)),
    d %>% filter(type == "binary") %>% transmute(ind, class, n, segment = "negative", pct = round(100 * neg / n, 1)),
    d %>% filter(type == "binary") %>% transmute(ind, class, n, segment = "positive", pct = round(100 * pos / n, 1))
  )
  long <- long %>% filter(!is.na(pct) & pct > 0)
  long$segment <- factor(long$segment, levels = c("negative", "low", "normal", "high", "positive"))
  long$ind_n <- paste0(long$ind, " (n=", long$n, ")")
  long$ind_n <- factor(long$ind_n, levels = rev(unique(long$ind_n)))
  lab <- long %>% group_by(ind_n) %>% slice_max(pct, n = 1, with_ties = FALSE) %>% ungroup()
  ggplot(long, aes(x = pct, y = ind_n, fill = segment)) +
    geom_col(width = bar_width, position = "stack") +
    #geom_text(data = lab, aes(label = pct), position = position_stack(vjust = 0.5), color = "white", size = 2.8) +
    geom_text(data = lab, aes(label = paste0(pct, "%")), position = position_stack(vjust = 0.5), color = "white", size = 2.8)+
    facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_fill_manual(values = colors, breaks = c("low", "high", "normal", "positive", "negative")) +
    labs(x = "Proportion (%)", y = NULL) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 7),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", size = 8, face = "bold"))
}
p_max <- plot_lab_stack(res, "Maximum disease activity")
p_max
p_mr <- plot_lab_stack(res, "Max response")
p_mr

######2 correlation analysis for clinical signatures and laboratory findings
#####(1) category variable correlation : auto-antibody (positive or negative)
library(openxlsx)
library(dplyr)
library(tidyr)

labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", colNames = FALSE)
ind <- as.character(unlist(labo[2, ]))
colnames(labo) <- trimws(sub("\\(.*", "", ind))
labo <- labo[-(1:2), ]
colnames(labo)[1:3] <- c("Number", "Name", "Stage")
colnames(labo)[colnames(labo) == "Serum amyloid A"] <- "SAA"
colnames(labo)[colnames(labo) == "Lymphocyte"] <- "CD4/CD8"
colnames(labo)[colnames(labo) == "Urine protein/creatinine ratio"] <- "Urine PCR"

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

mat <- lab_pick(labo, lab_ref)
mat_bin <- mat$binary

raw <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx", colNames = FALSE)
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
dim(clin_mat) 


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
    k <- k + 1
    out[[k]] <- data.frame(var1 = a, var2 = b, n = n,
                           rho = as.numeric(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE)
  }
  res <- bind_rows(out)
  res$p_adj <- p.adjust(res$p, method = "BH")
  res$dir <- ifelse(res$rho > 0, "↑", "↓")
  res[order(res$p_adj, res$p), ]
}

res_bin <- corr_binary(mat_bin, clin_mat)

library(psych)
library(pheatmap)
library(RColorBrewer)
colors2 <- colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
cor_heat <- function(x_mat, y_mat, sig_thr = 0.05, colors = colors2,
                     treeheight_row = 15, treeheight_col = 15, transpose = FALSE) {
  if (transpose) { tmp <- x_mat; x_mat <- y_mat; y_mat <- tmp }
  cor_res <- corr.test(x = x_mat, y = y_mat, method = "spearman", use = "pairwise", adjust = "none")
  r <- cor_res$r; p <- cor_res$p
  sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  keep_r <- rowSums(p < sig_thr, na.rm = TRUE) > 0
  keep_c <- colSums(p < sig_thr, na.rm = TRUE) > 0
  if (any(keep_r) & any(keep_c)) { r <- r[keep_r, keep_c, drop = FALSE]; p <- p[keep_r, keep_c, drop = FALSE]; sig <- sig[keep_r, keep_c, drop = FALSE] }
  sig[is.na(p)] <- ""
  max_abs <- max(abs(r), na.rm = TRUE); if (max_abs == 0) max_abs <- 1
  pheatmap(r, color = colorRampPalette(colors)(100),
           breaks = seq(-max_abs, max_abs, length.out = 100),
           cluster_rows = TRUE, cluster_cols = TRUE,
           treeheight_row = treeheight_row, treeheight_col = treeheight_col,
           display_numbers = sig, number_color = "black",
           fontsize = 7, border_color = NA, angle_col = 45)
  invisible(list(r = r, p = p))
}
p3 <- cor_heat(mat_bin, clin_mat)
p4 <- cor_heat(mat_bin, clin_mat, transpose = TRUE)

#####(2) numeric variable :  check higher or lower than reference range
library(openxlsx)
library(dplyr)
library(tidyr)
library(psych)
library(pheatmap)
library(RColorBrewer)

labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", colNames = FALSE)
ind <- as.character(unlist(labo[2, ]))
colnames(labo) <- trimws(sub("\\(.*", "", ind))
labo <- labo[-(1:2), ]
colnames(labo)[1:3] <- c("Number", "Name", "Stage")
colnames(labo)[colnames(labo) == "Serum amyloid A"] <- "SAA"
colnames(labo)[colnames(labo) == "Lymphocyte"] <- "CD4/CD8"
colnames(labo)[colnames(labo) == "Urine protein/creatinine ratio"] <- "Urine PCR"
labo$Stage <- trimws(as.character(labo$Stage))

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

stages <- c("Before genetic diagnosis", "Maximum disease activity", "Max response", "Last follow-up")
mat_num_up <- setNames(lapply(stages, function(st) lab_pick_binary(labo, lab_ref, st)), stages)
lapply(mat_num_up, dim)
lapply(mat_num_up, colnames)

raw <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Clinical data.xlsx", colNames = FALSE)
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
dim(clin_mat)

drop_clin <- c("Complement system dysregulation/activation",
               "SLE overlap / Systemic lupus erythematosus",
               "Positive autoantibodies","Crohn’s disease")
keep_clin <- !trimws(colnames(clin_mat)) %in% drop_clin
clin_mat <- clin_mat[, keep_clin, drop = FALSE]
dim(clin_mat)

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
  res$dir <- ifelse(res$rho > 0, "↑", "↓")
  res[order(res$p_adj, res$p), ]
}

res_num_up <- lapply(mat_num_up, function(mm) corr_binary(mm, clin_mat))
names(res_num_up) <- names(mat_num_up)
lapply(res_num_up, nrow)

colors2 <- colorRampPalette(c("#542788", "#F7F7F7", "#D6604D"))(101)
dist_pw <- function(M) {
  n <- nrow(M); d <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    ok <- !is.na(M[i, ]) & !is.na(M[j, ])
    d[i, j] <- if (sum(ok) > 0) sqrt(sum((M[i, ok] - M[j, ok])^2)) else NA
  }
  if (all(is.na(d))) d[] <- 0 else d[is.na(d)] <- max(d, na.rm = TRUE)
  as.dist(d)
}
cor_heat <- function(x_mat, y_mat, sig_thr = 0.05, colors = colors2,
                     treeheight_row = 10, treeheight_col = 10, transpose = FALSE) {
  if (transpose) { tmp <- x_mat; x_mat <- y_mat; y_mat <- tmp }
  cor_res <- corr.test(x = x_mat, y = y_mat, method = "spearman", use = "pairwise", adjust = "none")
  r <- cor_res$r; p <- cor_res$p
  sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  keep_r <- rowSums(p < sig_thr, na.rm = TRUE) > 0
  keep_c <- colSums(p < sig_thr, na.rm = TRUE) > 0
  if (any(keep_r) & any(keep_c)) { r <- r[keep_r, keep_c, drop = FALSE]; p <- p[keep_r, keep_c, drop = FALSE]; sig <- sig[keep_r, keep_c, drop = FALSE] }
  sig[is.na(p)] <- ""
  max_abs <- max(abs(r), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  pheatmap(r, color = colorRampPalette(colors)(100),
           breaks = seq(-max_abs, max_abs, length.out = 100),
           cluster_rows = nrow(r) > 1, cluster_cols = ncol(r) > 1,
           clustering_distance_rows = if (nrow(r) > 1) dist_pw(r) else NULL,
           clustering_distance_cols = if (ncol(r) > 1) dist_pw(t(r)) else NULL,
           treeheight_row = treeheight_row, treeheight_col = treeheight_col,
           display_numbers = sig, number_color = "black",
           fontsize = 7, border_color = NA, angle_col = 45)
  invisible(list(r = r, p = p))
}

cor_heat(mat_num_up[["Before genetic diagnosis"]], clin_mat)
cor_heat(mat_num_up[["Maximum disease activity"]], clin_mat)
cor_heat(mat_num_up[["Max response"]], clin_mat)
cor_heat(mat_num_up[["Last follow-up"]], clin_mat)

######3 Specific laboratory findings for 4 stages
plot_lab_panel <- function(labo, treat, pts, vars, line_col="#2E75B6", line_size=1, point_col="#C00000", point_size=2, stage_levels=c("Before genetic diagnosis","Maximum disease activity","Max response","Last follow-up"), stage_short=c("Before genetic diagnosis","Disease flare","Maximal response","Last follow-up")){
  labo <- as.data.frame(labo, stringsAsFactors = FALSE)
  if (!("Number" %in% names(labo)) && identical(trimws(as.character(labo[2,1])), "Number")) { names(labo) <- as.character(labo[2, ]); labo <- labo[-(1:2), ]; rownames(labo) <- NULL }
  for (v in c("Number","Stage")) labo[[v]] <- as.character(labo[[v]])
  norm <- function(x) gsub("[^[:alnum:]]","",tolower(x))
  nm <- names(labo)
  alias <- c(SAA="Serum amyloid A (SAA, <10 mg/L)","CD4/CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)","CD4.CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)",UPCR="Urine protein/creatinine ratio (0-30 mg/g)","Urine PCR"="Urine protein/creatinine ratio (0-30 mg/g)")
  idx <- match(vars, nm)
  if (any(is.na(idx))) { i2 <- match(norm(vars[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i2 }
  v2 <- ifelse(vars %in% names(alias), unname(alias[vars]), vars)
  if (any(is.na(idx))) { i3 <- match(v2[is.na(idx)], nm); idx[is.na(idx)] <- i3
  if (any(is.na(idx))) { i4 <- match(norm(v2[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i4 } }
  miss <- vars[is.na(idx)]
  if (length(miss)) stop("columns not found in labo: ", paste(miss, collapse=", "))
  use <- nm[idx]
  for (v in use) labo[[v]] <- suppressWarnings(as.numeric(as.character(labo[[v]])))
  labo <- labo[labo$Number %in% pts & labo$Stage %in% stage_levels, c("Number","Stage", use)]
  long <- pivot_longer(labo, all_of(use), names_to="Indicator", values_to="Value")
  long$Stage <- factor(long$Stage, levels=stage_levels)
  long$Indicator <- factor(long$Indicator, levels=use)
  dm <- data.frame(full=c("Prednisolone","Methylprednisolone","Methotrexate (MTX)","Mycophenolate mofetil (MMF)","Sulfasalazine","Cyclophosphamide","Leflunomide","Hydroxychloroquine (HCQ)","Cyclosporine","Tacrolimus","Adalimumab","Etanercept","Infliximab","Yisaipu","Anakinra","Canakinumab","Firsekibart","Tocilizumab","Belimumab","Tofacitinib","Baricitinib","Thalidomide","IVIG","Colchicine","Growth hormone"), abbr=c("PSL","mPSL","MTX","MMF","SSZ","CYC","LEF","HCQ","CsA","TAC","ADA","ETA","IFX","YSP","ANA","CAN","FIR","TCZ","BLM","TOF","BARI","THD","IVIG","COL","GH"), cat=c("GC","GC","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","Anti-TNF","Anti-TNF","Anti-TNF","Anti-TNF","Anti-IL-1","Anti-IL-1","Anti-IL-1","Anti-IL-6","B-cell","JAK","JAK","Immuno","IVIG","Colch","GH"), stringsAsFactors=FALSE)
  dm$key <- gsub("\\s+","",tolower(dm$full))
  cat_order <- c("GC","csDMARD","Anti-TNF","Anti-IL-1","Anti-IL-6","B-cell","JAK","Immuno","IVIG","Colch","GH")
  tt <- as.data.frame(treat, stringsAsFactors = FALSE)
  for (v in c("Number","Stage")) tt[[v]] <- as.character(tt[[v]])
  cn <- names(tt)
  pick <- function(pat) cn[grepl(pat, cn)][1]
  c_glu <- pick("^Glucocort"); c_dm <- pick("DMARD"); c_tnf <- pick("Anti.TNF"); c_il1 <- pick("Anti.IL.1"); c_il6 <- pick("Anti.IL.6"); c_b <- pick("B.cell"); c_jak <- pick("^JAK"); c_ivig <- pick("IVIG"); c_col <- pick("^Colch"); c_gh <- pick("Growth"); c_imm <- pick("^Immunomod"); c_oth <- pick("^Other.supplemental")
  parse_treat <- function(pid, stg) {
    r <- tt[tt$Number == pid & tt$Stage == stg, , drop=FALSE]
    if (!nrow(r)) return(NULL)
    toks <- character(0)
    take <- function(col, yesval=NULL) {
      if (is.na(col)) return(invisible(NULL))
      x <- trimws(as.character(r[[col]][1]))
      if (is.na(x) || x %in% c("", "unused", "unknown")) return(invisible(NULL))
      if (!is.null(yesval)) { if (x == "yes") toks <<- c(toks, yesval); return(invisible(NULL)) }
      if (x == "yes") return(invisible(NULL))
      for (pt in trimws(strsplit(x, ",")[[1]])) {
        k <- gsub("\\s+", "", tolower(pt))
        i <- match(k, dm$key)
        if (is.na(i)) i <- match(sub("\\(.*\\)$", "", k), dm$key)
        if (!is.na(i)) toks <<- c(toks, dm$abbr[i])
      }
      invisible(NULL)
    }
    take(c_glu); take(c_dm); take(c_tnf); take(c_il1); take(c_il6); take(c_b); take(c_jak); take(c_ivig, "IVIG"); take(c_col, "Colchicine"); take(c_gh, "Growth hormone"); take(c_imm); take(c_oth)
    if (!length(toks)) return(NULL)
    ab <- unique(toks)
    ord <- order(match(dm$cat[match(ab, dm$abbr)], cat_order))
    ab <- ab[ord]
    c(paste(ab, collapse="+"), paste(dm$full[match(ab, dm$abbr)], collapse="+"))
  }
  drug_rows <- lapply(pts, function(p) {
    m <- sapply(stage_levels, function(s) { r <- parse_treat(p, s); if (is.null(r)) c(NA_character_, NA_character_) else r })
    data.frame(Number=p, Stage=stage_levels, Drug_abbr=m[1, ], Drug_full=m[2, ], stringsAsFactors=FALSE)
  })
  drug_df <- bind_rows(drug_rows)
  drug_wide <- pivot_wider(drug_df, id_cols=Number, names_from=Stage, values_from=Drug_abbr)
  ind_lbl <- setNames(sub(" \\(", "\n(", use), use)
  plots <- list()
  for (p in pts) {
    dat <- long[long$Number == p, ]
    dd <- drug_df[drug_df$Number == p, ]
    xl <- sapply(seq_along(stage_levels), function(i) { d <- dd$Drug_abbr[i]; if (!is.na(d)) paste0(stage_short[i], "\n", d) else stage_short[i] })
    plots[[p]] <- ggplot(dat, aes(Stage, Value)) + geom_line(aes(group=1), color=line_col, linewidth=line_size) + geom_point(color=point_col, size=point_size) + facet_wrap(~ Indicator, ncol=1, scales="free_y", strip.position="left", labeller=labeller(Indicator=ind_lbl)) + scale_x_discrete(labels=xl) + labs(title=p) + theme_bw() + theme(axis.title=element_blank(), axis.text.x=element_text(size=5.5, lineheight=0.85), strip.text.y.left=element_text(size=7, angle=0, lineheight=0.9), strip.background=element_rect(fill="grey92", linewidth=0.2), plot.title=element_text(size=10, face="bold", hjust=0.5, margin=margin(b=2)), panel.spacing=unit(0.15, "lines"), panel.grid.minor=element_blank())
  }
  list(plots=plots, drugs=drug_wide, drug_long=drug_df)
}
labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", sheetIndex = 1, header = FALSE)
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", sheetIndex = 1)
pts <- c("P1","P3","P5","P7","P8","P9","P15","P16","P18","P20","P21","P25","P26");length(pts)
vars <- c("Serum creatinine (100-110 μmol/L)","AST (<40 U/L)","ALT (<50 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)","Serum amyloid A (SAA, <10 mg/L)")
res <- plot_lab_panel(labo, treat, pts=pts, vars=vars, line_col="#2E75B6", line_size=1.2, point_col="#C00000", point_size=2.5)
names(res$plots)
res$drugs
res$plots$P16
length(res$drugs)

library(ggpubr)
library(gridExtra)
n_p <- length(res$plots)
setwd("E:\\Cohort PPT\\JIA\\code")
pdf("lab_panel_all.pdf", width = 13.5, height = 9.5, onefile = TRUE)
for (i in seq_len(ceiling(n_p / 4))) {
  idx <- ((i - 1) * 4 + 1):min(i * 4, n_p)
  grid.arrange(grobs = res$plots[idx], ncol = 4, respect = TRUE)
}
dev.off()


suppressMessages({library(xlsx);library(dplyr);library(tidyr);library(ggplot2)})
plot_lab_panel <- function(labo, treat, pts, vars, line_col="#2E75B6", line_size=1, point_size=2, col_high="red", col_low="blue", col_norm="black", stage_levels=c("Before genetic diagnosis","Maximum disease activity","Max response","Last follow-up"), stage_short=c("Before genetic diagnosis","Disease flare","Maximal response","Last follow-up")){
  labo <- as.data.frame(labo, stringsAsFactors = FALSE)
  if (!("Number" %in% names(labo)) && identical(trimws(as.character(labo[2,1])), "Number")) { names(labo) <- as.character(labo[2, ]); labo <- labo[-(1:2), ]; rownames(labo) <- NULL }
  for (v in c("Number","Stage")) labo[[v]] <- as.character(labo[[v]])
  norm <- function(x) gsub("[^[:alnum:]]","",tolower(x))
  nm <- names(labo)
  alias <- c(SAA="Serum amyloid A (SAA, <10 mg/L)","CD4/CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)","CD4.CD8"="Lymphocyte (CD4+T / CD8+T ,0.9-2.0)",UPCR="Urine protein/creatinine ratio (0-30 mg/g)","Urine PCR"="Urine protein/creatinine ratio (0-30 mg/g)")
  idx <- match(vars, nm)
  if (any(is.na(idx))) { i2 <- match(norm(vars[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i2 }
  v2 <- ifelse(vars %in% names(alias), unname(alias[vars]), vars)
  if (any(is.na(idx))) { i3 <- match(v2[is.na(idx)], nm); idx[is.na(idx)] <- i3
  if (any(is.na(idx))) { i4 <- match(norm(v2[is.na(idx)]), norm(nm)); idx[is.na(idx)] <- i4 } }
  miss <- vars[is.na(idx)]
  if (length(miss)) stop("columns not found in labo: ", paste(miss, collapse=", "))
  use <- nm[idx]
  for (v in use) labo[[v]] <- suppressWarnings(as.numeric(as.character(labo[[v]])))
  labo <- labo[labo$Number %in% pts & labo$Stage %in% stage_levels, c("Number","Stage", use)]
  long <- pivot_longer(labo, all_of(use), names_to="Indicator", values_to="Value")
  long$Stage <- factor(long$Stage, levels=stage_levels)
  long$Indicator <- factor(long$Indicator, levels=use)
  ref <- data.frame(Indicator=c("Serum creatinine (100-110 μmol/L)","AST (<40 U/L)","ALT (<50 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)","Serum amyloid A (SAA, <10 mg/L)"), lo=c(100,0,0,0.068,0,0), hi=c(110,40,50,8.2,20,10), stringsAsFactors=FALSE)
  long <- left_join(long, ref, by="Indicator")
  long$lo[is.na(long$lo)] <- -Inf; long$hi[is.na(long$hi)] <- Inf
  long$Status <- ifelse(long$Value > long$hi, "high", ifelse(long$Value < long$lo, "low", "normal"))
  dm <- data.frame(full=c("Prednisolone","Methylprednisolone","Methotrexate (MTX)","Mycophenolate mofetil (MMF)","Sulfasalazine","Cyclophosphamide","Leflunomide","Hydroxychloroquine (HCQ)","Cyclosporine","Tacrolimus","Adalimumab","Etanercept","Infliximab","Yisaipu","Anakinra","Canakinumab","Firsekibart","Tocilizumab","Belimumab","Tofacitinib","Baricitinib","Thalidomide","IVIG","Colchicine","Growth hormone"), abbr=c("PSL","mPSL","MTX","MMF","SSZ","CYC","LEF","HCQ","CsA","TAC","ADA","ETA","IFX","YSP","ANA","CAN","FIR","TCZ","BLM","TOF","BARI","THD","IVIG","COL","GH"), cat=c("GC","GC","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","csDMARD","Anti-TNF","Anti-TNF","Anti-TNF","Anti-TNF","Anti-IL-1","Anti-IL-1","Anti-IL-1","Anti-IL-6","B-cell","JAK","JAK","Immuno","IVIG","Colch","GH"), stringsAsFactors=FALSE)
  dm$key <- gsub("\\s+","",tolower(dm$full))
  cat_order <- c("GC","csDMARD","Anti-TNF","Anti-IL-1","Anti-IL-6","B-cell","JAK","Immuno","IVIG","Colch","GH")
  tt <- as.data.frame(treat, stringsAsFactors = FALSE)
  for (v in c("Number","Stage")) tt[[v]] <- as.character(tt[[v]])
  cn <- names(tt)
  pick <- function(pat) cn[grepl(pat, cn)][1]
  c_glu <- pick("^Glucocort"); c_dm <- pick("DMARD"); c_tnf <- pick("Anti.TNF"); c_il1 <- pick("Anti.IL.1"); c_il6 <- pick("Anti.IL.6"); c_b <- pick("B.cell"); c_jak <- pick("^JAK"); c_ivig <- pick("IVIG"); c_col <- pick("^Colch"); c_gh <- pick("Growth"); c_imm <- pick("^Immunomod"); c_oth <- pick("^Other.supplemental")
  parse_treat <- function(pid, stg) {
    r <- tt[tt$Number == pid & tt$Stage == stg, , drop=FALSE]
    if (!nrow(r)) return(NULL)
    toks <- character(0)
    take <- function(col, yesval=NULL) {
      if (is.na(col)) return(invisible(NULL))
      x <- trimws(as.character(r[[col]][1]))
      if (is.na(x) || x %in% c("", "unused", "unknown")) return(invisible(NULL))
      if (!is.null(yesval)) { if (x == "yes") toks <<- c(toks, yesval); return(invisible(NULL)) }
      if (x == "yes") return(invisible(NULL))
      for (pt in trimws(strsplit(x, ",")[[1]])) {
        k <- gsub("\\s+", "", tolower(pt))
        i <- match(k, dm$key)
        if (is.na(i)) i <- match(sub("\\(.*\\)$", "", k), dm$key)
        if (!is.na(i)) toks <<- c(toks, dm$abbr[i])
      }
      invisible(NULL)
    }
    take(c_glu); take(c_dm); take(c_tnf); take(c_il1); take(c_il6); take(c_b); take(c_jak); take(c_ivig, "IVIG"); take(c_col, "Colchicine"); take(c_gh, "Growth hormone"); take(c_imm); take(c_oth)
    if (!length(toks)) return(NULL)
    ab <- unique(toks)
    ord <- order(match(dm$cat[match(ab, dm$abbr)], cat_order))
    ab <- ab[ord]
    c(paste(ab, collapse="+"), paste(dm$full[match(ab, dm$abbr)], collapse="+"))
  }
  drug_rows <- lapply(pts, function(p) {
    m <- sapply(stage_levels, function(s) { r <- parse_treat(p, s); if (is.null(r)) c(NA_character_, NA_character_) else r })
    data.frame(Number=p, Stage=stage_levels, Drug_abbr=m[1, ], Drug_full=m[2, ], stringsAsFactors=FALSE)
  })
  drug_df <- bind_rows(drug_rows)
  drug_wide <- pivot_wider(drug_df, id_cols=Number, names_from=Stage, values_from=Drug_abbr)
  ind_lbl <- setNames(sub(" \\(", "\n(", use), use)
  plots <- list()
  for (p in pts) {
    dat <- long[long$Number == p, ]
    dd <- drug_df[drug_df$Number == p, ]
    xl <- sapply(seq_along(stage_levels), function(i) { d <- dd$Drug_abbr[i]; if (!is.na(d)) paste0(stage_short[i], "\n", d) else stage_short[i] })
    plots[[p]] <- ggplot(dat, aes(Stage, Value)) + geom_line(aes(group=1), color=line_col, linewidth=line_size) + geom_point(aes(color=Status), size=point_size) + facet_wrap(~ Indicator, ncol=1, scales="free_y", strip.position="top", labeller=labeller(Indicator=ind_lbl)) + scale_color_manual(values=c(high=col_high, normal=col_norm, low=col_low), name=NULL, breaks=c("high","normal","low"), labels=c("Above normal","Normal","Below normal")) + scale_x_discrete(labels=xl) + labs(title=p) + theme_classic() + theme(axis.title=element_blank(), axis.text.x=element_text(size=5.5, lineheight=0.85), axis.text.y=element_text(size=6), axis.line=element_line(color="grey40", linewidth=0.3), panel.grid.major.x=element_line(linetype="22", color="grey85", linewidth=0.4), panel.grid.major.y=element_line(linetype="11", color="grey85", linewidth=0.4), panel.grid.minor=element_blank(), strip.background=element_blank(), strip.text=element_text(face="bold", size=9), plot.title=element_text(size=11, face="bold", hjust=0.5, margin=margin(b=2)), panel.spacing=unit(0.2, "lines"), legend.position="top", legend.title=element_blank(), legend.text=element_text(size=7))
  }
  list(plots=plots, drugs=drug_wide, drug_long=drug_df)
}
labo <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\Labo.xlsx", sheetIndex = 1, header = FALSE)
treat <- read.xlsx("E:\\Cohort PPT\\JIA\\code\\treatment data.xlsx", sheetIndex = 1)
pts <- c("P5","P21")
vars <- c("AST (<40 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)")

names(res$plots)

res <- plot_lab_panel(labo, treat, pts = c("P5","P21","P25"),
                      vars = c("AST (<40 U/L)","CRP (0.068-8.2 mg/L)","ESR (0-20mm/h)"),
                      line_col = "gray60", line_size = 1.2,
                      col_high = "#BB0021", 
                      col_norm = "black", 
                      col_low = "#6699CC",
                      point_size = 3.5)
library(patchwork)
res$plots$P5|res$plots$P21


