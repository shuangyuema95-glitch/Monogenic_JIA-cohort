#####GSE80060, sJIA whole blood, canakinumab, GPL570
dir <- "E:/Cohort PPT/JIA/bulkdata"; setwd(dir)

########1. preprocessing
library(limma)
library(readxl)

prep_gse <- function(fmat, fplat, transform = c("none", "log2"), qn = TRUE, cut = 5, frac = 0.2,
                     keep_table = NULL) {
  transform <- match.arg(transform)
  
  rd_expr <- function(f) {
    ln <- readLines(f)
    tb <- ln[(grep("^!series_matrix_table_begin", ln) + 1):(grep("^!series_matrix_table_end", ln) - 1)]
    con <- textConnection(tb)
    e <- as.matrix(read.delim(con, row.names = 1, check.names = FALSE))
    close(con)
    if (transform == "log2" && max(e, na.rm = TRUE) > 20) e <- log2(e + 1)
    if (qn) e <- normalizeBetweenArrays(e, method = "quantile")
    e
  }
  
  rd_plat <- function(f) {
    pl <- readLines(f)
    pb <- grep("^!platform_table_begin", pl)
    pe <- grep("^!platform_table_end", pl)
    pt <- if (length(pb)) pl[(pb + 1):(pe - 1)] else pl[grep("^ID\t", pl)[1]:length(pl)]
    ph <- strsplit(pt[1], "\t", fixed = TRUE)[[1]]
    sx <- grep("gene symbol", ph, ignore.case = TRUE)[1]
    x2 <- strsplit(pt[-1], "\t", fixed = TRUE)
    s2 <- lapply(x2, function(z) {
      zz <- if (is.na(z[sx])) "" else z[sx]
      trimws(strsplit(zz, "///", fixed = TRUE)[[1]])
    })
    pid <- vapply(x2, `[`, "", 1)
    data.frame(PROBEID = rep(pid, lengths(s2)), SYMBOL = trimws(unlist(s2)), stringsAsFactors = FALSE)
  }
  
  e <- rd_expr(fmat[1])
  m <- rd_plat(fplat[1])
  info <- keep_table
  info$col <- info$gsm
  info <- info[info$col %in% colnames(e), ]
  info <- info[!duplicated(info$col), ]
  e <- e[, info$col, drop = FALSE]
  
  sel <- !is.na(m$SYMBOL) & m$SYMBOL != "" & m$SYMBOL != "---" & m$PROBEID %in% rownames(e)
  m <- m[sel, ]
  m <- m[order(-rowMeans(e)[m$PROBEID]), ]
  m <- m[!duplicated(m$SYMBOL), ]
  
  e2 <- e[m$PROBEID, , drop = FALSE]
  rownames(e2) <- m$SYMBOL
  e2 <- e2[rowSums(e2 > cut) >= ceiling(frac * ncol(e2)), , drop = FALSE]
  list(expr = e2, info = info)
}

raw <- as.data.frame(read_excel("GSE80060_JIAcana.xlsx", sheet = "Sheet2", col_names = FALSE))
ch <- as.data.frame(lapply(raw[, 3:8], function(x) sub("^[^:]+:\\s*", "", as.character(x))))
names(ch) <- c("tissue", "disease", "subject", "visit", "treatment", "acr")

sel <- data.frame(gsm = as.character(raw[[2]]),
                  title = as.character(raw[[1]]),
                  ch,
                  stringsAsFactors = FALSE)
sel <- sel[grepl("^GSM", sel$gsm), ]
sel$acr <- suppressWarnings(as.numeric(sel$acr))
sel <- sel[sel$disease != "SJIA" | sel$treatment == "Canakinumab", ]

g2 <- prep_gse("GSE80060_series_matrix.txt.gz", "GPL570-55999.txt",
               transform = "log2", qn = TRUE, cut = 5, frac = 0.2,
               keep_table = sel)

expr <- g2$expr
info2 <- g2$info

info2$study <- ifelse(grepl("^SJIA_[0-9]_", info2$subject),
                      sub("^SJIA_([0-9])_.*$", "\\1", info2$subject), NA)
info2$group <- factor(ifelse(info2$disease == "SJIA", "sJIA", "HC"), c("HC", "sJIA"))
info2$grp3 <- factor(ifelse(info2$group == "HC", "HC",
                            ifelse(info2$visit == "Day1_BL", "BL", "D3")), c("HC", "BL", "D3"))
info2$resp <- ifelse(is.na(info2$acr), NA, ifelse(info2$acr >= 70, "R", "NR"))

info2 <- info2[order(info2$grp3, info2$subject, info2$visit), ]
expr <- expr[, info2$col, drop = FALSE]

dim(expr)
range(expr)
sum(is.na(expr))
table(info2$grp3)
table(info2$grp3, info2$resp, useNA = "ifany")


########2. cell infiltration
library(ggplot2)
####(1) deconvolution
lin <- 2^expr - 1
write.table(round(lin, 3), "GSE80060_mix.txt", sep = "\t", quote = FALSE, col.names = NA)
source("CIBERSORT.R")
res <- CIBERSORT("LM22.txt", "GSE80060_mix.txt", perm = 100, QN = FALSE)

f <- res[, 1:22, drop = FALSE]
info2$grpB <- factor(ifelse(info2$grp3 == "HC", "HC",
                            ifelse(info2$grp3 == "D3", info2$resp, NA)), c("HC", "R", "NR"))
table(res[, "P-value"] < 0.05)
table(info2$grpB, useNA = "ifany")

####(2) pairwise tests per cell type
cel_test <- function(f, grp, lev, lab) {
  d <- data.frame(cell = rep(colnames(f), each = nrow(f)),
                  val = as.vector(as.matrix(f)),
                  grp = rep(as.character(grp), times = ncol(f)),
                  stringsAsFactors = FALSE)
  d$cell <- factor(d$cell, levels = colnames(f))
  d$grp <- factor(d$grp, levels = lev, labels = lab)
  d <- d[!is.na(d$grp), ]
  cb <- list(c(1, 2), c(1, 3), c(2, 3))
  out <- lapply(levels(d$cell), function(k) {
    s <- d[d$cell == k, ]
    do.call(rbind, lapply(seq_along(cb), function(j) {
      a <- cb[[j]][1]; b <- cb[[j]][2]
      s2 <- s[as.integer(s$grp) %in% c(a, b), ]
      s2$grp <- droplevels(s2$grp)
      p <- if (max(s$val) == 0 || nlevels(s2$grp) < 2) NA else
        tryCatch(wilcox.test(val ~ grp, data = s2)$p.value, error = function(e) NA)
      data.frame(cell = k, pair = paste(lab[a], "vs", lab[b]),
                 n1 = sum(s2$grp == lab[a]), n2 = sum(s2$grp == lab[b]),
                 m1 = mean(s2$val[s2$grp == lab[a]]), m2 = mean(s2$val[s2$grp == lab[b]]),
                 x = a, xend = b, p = p, stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  out$padj <- ave(out$p, out$pair, FUN = function(z) p.adjust(z, "BH"))
  out$star <- as.character(cut(out$padj, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "")))
  out$star[is.na(out$star)] <- "ns"
  out
}

####(3) plots for selected cell types
pl_cell <- function(f, grp, cols, lev, lab, cells, tab,
                    pt_size = 2.1, box_w = 0.6, fill_alpha = 0.85, pair = NULL) {
  dk <- function(x, k = 0.62) {
    m <- col2rgb(x) / 255
    rgb(m[1] * k, m[2] * k, m[3] * k)
  }
  nn <- table(factor(as.character(grp), levels = lev, labels = lab))
  xl <- setNames(paste0(lab, "\n(n=", as.integer(nn), ")"), lab)
  d <- data.frame(cell = rep(colnames(f), each = nrow(f)),
                  val = as.vector(as.matrix(f)),
                  grp = rep(as.character(grp), times = ncol(f)),
                  stringsAsFactors = FALSE)
  d$cell <- factor(d$cell, levels = colnames(f))
  d$grp <- factor(d$grp, levels = lev, labels = lab)
  d <- d[!is.na(d$grp), ]
  dc <- sapply(cols, dk)
  if (is.null(cells)) {
    t3 <- tab[tab$padj < 0.05 & !is.na(tab$padj), ]
    if (!is.null(pair)) t3 <- t3[t3$pair %in% pair, ]
    cells <- as.character(unique(t3$cell))
  }
  out <- lapply(setNames(cells, cells), function(k) {
    s <- d[d$cell == k, ]
    t2 <- tab[tab$cell == k & tab$star != "" & !is.na(tab$star), ]
    if (!is.null(pair)) t2 <- t2[t2$pair %in% pair, ]
    g <- ggplot(s, aes(grp, val)) +
      geom_boxplot(aes(fill = grp), width = box_w, outlier.shape = NA, linewidth = 0.3, colour = "black") +
      geom_point(aes(colour = grp), shape = 19, size = pt_size, alpha = 0.85,
                 position = position_jitterdodge(jitter.width = box_w * 0.2, dodge.width = box_w)) +
      scale_fill_manual(values = ggplot2::alpha(cols, fill_alpha)) +
      scale_colour_manual(values = dc) +
      scale_x_discrete(labels = xl) +
      guides(colour = "none", fill = "none") +
      labs(x = NULL, y = "Estimated proportion", title = k) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 10, colour = "black", hjust = 0.5),
            legend.position = "none")
    if (nrow(t2)) {
      t2 <- t2[order(t2$x), ]
      y0 <- max(s$val)
      t2$y <- y0 * (1.05 + 0.09 * (seq_len(nrow(t2)) - 1))
      t2$ty <- t2$y * 1.03
      g <- g +
        geom_segment(data = t2, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
        geom_text(data = t2, aes(x = (x + xend) / 2, y = ty, label = star), inherit.aes = FALSE, size = 3, vjust = 0) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.3)))
    } else {
      g <- g + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    }
    g
  })
  out
}

####(4) comparisons and plots
tabA <- cel_test(f, info2$grp3, c("HC", "BL", "D3"), c("HC", "BL", "D3"))
tabB <- cel_test(f, info2$grpB, c("HC", "R", "NR"), c("HC", "R", "NR"))

tabA[tabA$padj < 0.05, ]
tabB[tabB$padj < 0.05, ]

plA <- pl_cell(f, info2$grp3, c(HC = "grey90", BL = "#B5BBE3", D3 = "#FEC260"),
               c("HC", "BL", "D3"), c("HC", "BL", "D3"), NULL, tabA)
plB <- pl_cell(f, info2$grpB, c(HC = "grey90", R = "#6699CC", NR = "#DF5734"),
               c("HC", "R", "NR"), c("HC", "R", "NR"), NULL, tabB)
names(plA)
names(plB)

library(patchwork)
pA_all <- if (length(plA)) wrap_plots(plA, ncol = min(4, length(plA))) else NULL
pB_all <- if (length(plB)) wrap_plots(plB, ncol = min(4, length(plB))) else NULL

pA_all
pB_all

plA$`B cells naive`#using this
plA$`T cells CD8`#using this
save(tabA,tabB,info2,file="E:/Cohort PPT/JIA/bulkdata/GSE80060_cellinfiltra.Rdata")

########3. differential expression
library(limma)
de <- list()

####(1) BL vs HC
sb <- info2[info2$grp3 %in% c("HC", "BL"), ]
sb$grp3 <- droplevels(sb$grp3)
design <- model.matrix(~ 0 + grp3, sb)
colnames(design) <- levels(sb$grp3)
cm <- makeContrasts(BL_HC = BL - HC, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr[, sb$gsm], design), cm))
de$BL_HC <- topTable(fit, coef = "BL_HC", number = Inf, sort.by = "P")

####(2) D3 vs HC
sd2 <- info2[info2$grp3 %in% c("HC", "D3"), ]
sd2$grp3 <- droplevels(sd2$grp3)
design <- model.matrix(~ 0 + grp3, sd2)
colnames(design) <- levels(sd2$grp3)
cm <- makeContrasts(D3_HC = D3 - HC, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr[, sd2$gsm], design), cm))
de$D3_HC <- topTable(fit, coef = "D3_HC", number = Inf, sort.by = "P")

####(3) D3 vs BL, paired
pp <- intersect(info2$subject[info2$visit == "Day1_BL" & info2$treatment == "Canakinumab"],
                info2$subject[info2$visit == "Day3" & info2$treatment == "Canakinumab"])
sp <- info2[info2$subject %in% pp & info2$treatment == "Canakinumab", ]
sp$subject <- factor(sp$subject)
sp$visit <- factor(sp$visit, c("Day1_BL", "Day3"))
design <- model.matrix(~ subject + visit, sp)
fit <- eBayes(lmFit(expr[, sp$gsm], design))
de$BL_D3 <- topTable(fit, coef = "visitDay3", number = Inf, sort.by = "P")

####(4) R vs NR at D3
sr <- info2[info2$grp3 == "D3" & !is.na(info2$resp), ]
sr$resp <- factor(sr$resp, c("R", "NR"))
sr$study <- factor(sr$study, c("1", "2"))
design <- model.matrix(~ 0 + resp + study, sr)
cm <- makeContrasts(R_NR = respR - respNR, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr[, sr$gsm], design), cm))
de$R_NR_D3 <- topTable(fit, coef = "R_NR", number = Inf, sort.by = "P")

####(5) direction check
chk <- c("ISG15", "IFI44L", "MX1", "S100A8", "S100A9", "FCGR3B", "CXCR2")
gm <- function(g, lev) sapply(g, function(x) mean(expr[x, info2$grp3 == lev]))

data.frame(gene = chk,
           t_BL_HC = de$BL_HC$t[match(chk, rownames(de$BL_HC))],
           d_BL_HC = gm(chk, "BL") - gm(chk, "HC"),
           t_BL_D3 = de$BL_D3$t[match(chk, rownames(de$BL_D3))],
           d_BL_D3 = gm(chk, "D3") - gm(chk, "BL"))

####(6) summary
smry <- do.call(rbind, lapply(names(de), function(k) {
  d <- de[[k]]
  s <- d[d$adj.P.Val < 0.05 & abs(d$logFC) > 1, ]
  data.frame(contrast = k,
             n_up = sum(s$logFC > 1), n_dn = sum(s$logFC < -1),
             n_p05 = sum(d$P.Value < 0.05),
             IFN28 = sum(rownames(s) %in% IFN28),
             IFN28_up = sum(rownames(s)[s$logFC > 0] %in% IFN28),
             NFKB = sum(rownames(s) %in% NFKB$gene),
             NFKB_up = sum(rownames(s)[s$logFC > 0] %in% NFKB$gene),
             stringsAsFactors = FALSE)
}))
smry




########4. differential expression
library(limma)

de <- list()

####(1) BL vs HC
sb <- info2[info2$grp3 %in% c("HC", "BL"), ]
sb$grp3 <- droplevels(sb$grp3)
design <- model.matrix(~ 0 + grp3, sb)
colnames(design) <- levels(sb$grp3)
cm <- makeContrasts(BL_HC = BL - HC, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr[, sb$gsm], design), cm))
de$BL_HC <- topTable(fit, coef = "BL_HC", number = Inf, sort.by = "P")

####(2) D3 vs HC
sd2 <- info2[info2$grp3 %in% c("HC", "D3"), ]
sd2$grp3 <- droplevels(sd2$grp3)
design <- model.matrix(~ 0 + grp3, sd2)
colnames(design) <- levels(sd2$grp3)
cm <- makeContrasts(D3_HC = D3 - HC, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr[, sd2$gsm], design), cm))
de$D3_HC <- topTable(fit, coef = "D3_HC", number = Inf, sort.by = "P")

####(3) D3 vs BL, paired
pp <- intersect(info2$subject[info2$visit == "Day1_BL" & info2$treatment == "Canakinumab"],
                info2$subject[info2$visit == "Day3" & info2$treatment == "Canakinumab"])
sp <- info2[info2$subject %in% pp & info2$treatment == "Canakinumab", ]
sp$subject <- factor(sp$subject)
sp$visit <- factor(sp$visit, c("Day1_BL", "Day3"))
design <- model.matrix(~ subject + visit, sp)
fit <- eBayes(lmFit(expr[, sp$gsm], design))
de$BL_D3 <- topTable(fit, coef = "visitDay3", number = Inf, sort.by = "P")

####(4) R vs NR at D3
sr <- info2[info2$grp3 == "D3" & !is.na(info2$resp), ]
sr$resp <- factor(sr$resp, c("R", "NR"))
sr$study <- factor(sr$study, c("1", "2"))
design <- model.matrix(~ 0 + resp + study, sr)
cm <- makeContrasts(R_NR = respR - respNR, levels = design)
fit <- eBayes(contrasts.fit(lmFit(expr[, sr$gsm], design), cm))
de$R_NR_D3 <- topTable(fit, coef = "R_NR", number = Inf, sort.by = "P")

####(5) direction check
chk <- c("ISG15", "IFI44L", "MX1", "S100A8", "S100A9", "FCGR3B", "CXCR2")
gm <- function(g, lev) sapply(g, function(x) mean(expr[x, info2$grp3 == lev]))

data.frame(gene = chk,
           t_BL_HC = de$BL_HC$t[match(chk, rownames(de$BL_HC))],
           d_BL_HC = gm(chk, "BL") - gm(chk, "HC"),
           t_BL_D3 = de$BL_D3$t[match(chk, rownames(de$BL_D3))],
           d_BL_D3 = gm(chk, "D3") - gm(chk, "BL"))

####(6) summary
smry <- do.call(rbind, lapply(names(de), function(k) {
  d <- de[[k]]
  s <- d[d$adj.P.Val < 0.05 & abs(d$logFC) > 1, ]
  data.frame(contrast = k,
             n_up = sum(s$logFC > 1), n_dn = sum(s$logFC < -1),
             n_p05 = sum(d$P.Value < 0.05),
             IFN28 = sum(rownames(s) %in% IFN28),
             IFN28_up = sum(rownames(s)[s$logFC > 0] %in% IFN28),
             NFKB = sum(rownames(s) %in% NFKB$gene),
             NFKB_up = sum(rownames(s)[s$logFC > 0] %in% NFKB$gene),
             stringsAsFactors = FALSE)
}))
smry


########4. GSEA
library(msigdbr)
library(clusterProfiler)

h <- msigdbr(species = "Homo sapiens", category = "H")
k <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
sets <- rbind(data.frame(term = sub("^HALLMARK_", "", h$gs_name), gene = h$gene_symbol),
              data.frame(term = sub("^KEGG_", "", k$gs_name), gene = k$gene_symbol))

gse <- lapply(de, function(x) {
  v <- setNames(x[["t"]], rownames(x))
  GSEA(sort(v[!is.na(v)], decreasing = TRUE), TERM2GENE = sets,
       pvalueCutoff = 0.2, verbose = FALSE)
})

lab <- c(BL_HC = "higher_in_BL", D3_HC = "higher_in_D3",
         BL_D3 = "higher_at_D3", R_NR_D3 = "higher_in_R")

gsd <- function(x, nm) {
  if (is.null(x)) return(NULL)
  d <- as.data.frame(x)
  d <- d[d$p.adjust < 0.05, c("ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment")]
  if (!nrow(d)) return(NULL)
  d$contrast <- nm
  d$dir <- ifelse(d$NES > 0, lab[nm], sub("higher_", "lower_", lab[nm]))
  d[order(-d$NES), ]
}

gdp <- do.call(rbind, lapply(names(gse), function(n) gsd(gse[[n]], n)))
gdp
dim(gdp)
write.xlsx(gdp,file="E:/Cohort PPT/JIA/bulkdata/GSE80060GSEAresult.xlsx")

########6. GSEA plots
library(enrichplot)
ann <- c(BL_D3 = "Canakinumab Day3 vs BL",
         R_NR_D3 = "Responders vs non-responders, Day3")
colr <- c(BL_D3 = "#4466B0", R_NR_D3 = "#DF5734")

sel <- data.frame(contrast = c("BL_D3", "BL_D3", "R_NR_D3", "R_NR_D3", "R_NR_D3", "R_NR_D3"),
                  id = c("INFLAMMATORY_RESPONSE", "TNFA_SIGNALING_VIA_NFKB",
                         "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE",
                         "INTERFERON_GAMMA_RESPONSE", "INTERFERON_ALPHA_RESPONSE"),
                  stringsAsFactors = FALSE)

gp <- function(cid, id, bs = 8) {
  r <- gdp[gdp$contrast == cid & gdp$ID == id, ]
  if (!nrow(r)) return(NULL)
  r <- r[which.min(r$p.adjust), ]
  ttl <- sprintf("%s\n%s   NES = %.2f, FDR = %.2g", id, ann[cid], r$NES, r$p.adjust)
  gseaplot2(gse[[cid]], geneSetID = id, title = ttl,
            color = unname(colr[cid]), pvalue_table = FALSE, base_size = bs,
            subplots = 1:2, rel_heights = c(1.5, 0.5))}
gl <- lapply(seq_len(nrow(sel)), function(i) gp(sel$contrast[i], sel$id[i]))
names(gl) <- paste(sel$contrast, sel$id, sep = "_")
gl


#####################
#####################merge results
##### 1 merge GSEA results of GSE23687
########GSE23687 GSEA panels, from saved result
library(enrichplot)
library(patchwork)

gse_m <- readRDS("E:/Cohort PPT/JIA/bulkdata/1112_gse.RDS")

lab_m <- c(pre_vs_post_all11 = "higher_after_MTX", R_vs_NR_post = "higher_in_R")

gsd_m <- function(x, nm) {
  if (is.null(x)) return(NULL)
  d <- as.data.frame(x)
  d <- d[d$p.adjust < 0.05, c("ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment")]
  if (!nrow(d)) return(NULL)
  d$contrast <- nm
  d$dir <- ifelse(d$NES > 0, lab_m[nm], sub("higher_", "lower_", lab_m[nm]))
  d[order(-d$NES), ]
}
gdp_m <- do.call(rbind, lapply(names(gse_m), function(n) gsd_m(gse_m[[n]], n)))

ann_m <- c(pre_vs_post_all11 = "MTX post vs pre, paired")
colr_m <- c(pre_vs_post_all11 = "#4466B0")

gpm <- function(id, cid = "pre_vs_post_all11", bs = 8) {
  r <- gdp_m[gdp_m$contrast == cid & gdp_m$ID == id, ]
  if (!nrow(r)) return(NULL)
  r <- r[which.min(r$p.adjust), ]
  ttl <- sprintf("%s\n%s   NES = %.2f, FDR = %.2g", id, ann_m[cid], r$NES, r$p.adjust)
  gseaplot2(gse_m[[cid]], geneSetID = id, title = ttl, color = colr_m[cid],
            pvalue_table = FALSE, base_size = bs, subplots = 1:2, rel_heights = c(1.5, 0.5))
}

gl_mtx <- list(TNFA_SIGNALING_VIA_NFKB = gpm("TNFA_SIGNALING_VIA_NFKB"),
               INFLAMMATORY_RESPONSE = gpm("INFLAMMATORY_RESPONSE"))
gl_mtx[[1]]# NES <0 ,lower in post
gl_mtx[[2]]# NES <0 ,lower in post

library(patchwork)
library(grid)

gl8 <- gl8[!vapply(gl8, is.null, TRUE)]
gl8g <- lapply(gl8, function(p) grid::grid.grabExpr(print(p), wrap.grobs = TRUE))

p8 <- wrap_plots(gl8g, ncol = 4,byrow = FALSE)
p8
saveRDS(p8,file="E:/Cohort PPT/JIA/bulkdata/GSEA result merge.RDS")

##### 2 merge cellular infiltration
library(cowplot)#plA form GSE80060, plGSE8650
save(plA,plGSE8650,file="E:/Cohort PPT/JIA/bulkdata/GSE80060_8650Cellinfil.Rdata")
plot_grid(plA$`B cells naive`,plA$`T cells CD8`
          ,plGSE8650$`T cells CD4 naive`,plGSE8650$`T cells CD4 memory resting`,
          ncol=4,align="hv")

