library(limma)
library(openxlsx)
library(ggplot2)
setwd("E:/Cohort PPT/JIA/bulkdata")

######1 Preprocessing
library(limma)
library(readxl)

prep_gse <- function(fmat, fplat, transform = c("none", "log2"), qn = TRUE, cut = 5, frac = 0.2,
                     keep_field = NULL, keep_value = NULL, keep_table = NULL) {
  if (length(fmat) != length(fplat)) stop("fmat and fplat length mismatch")
  transform <- match.arg(transform)
  
  rd_expr <- function(f) {
    ln <- readLines(f)
    tb <- ln[(grep("^!series_matrix_table_begin", ln) + 1):(grep("^!series_matrix_table_end", ln) - 1)]
    con <- textConnection(tb)
    e <- as.matrix(read.delim(con, row.names = 1, check.names = FALSE))
    close(con)
    if (transform == "log2") e <- log2(e + 1)
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
  
  e <- lapply(fmat, rd_expr)
  pmap <- NULL
  if (length(e) == 2) {
    ia <- as.integer(sub("GSM", "", colnames(e[[1]])))
    ib <- as.integer(sub("GSM", "", colnames(e[[2]])))
    j <- match(ib - 1, ia)
    if (any(is.na(j))) stop("GSM pairing failed")
    pmap <- setNames(colnames(e[[1]])[j], colnames(e[[2]]))
    colnames(e[[2]]) <- unname(pmap)
    e <- rbind(e[[1]], e[[2]])
  } else {
    e <- e[[1]]
  }
  m <- do.call(rbind, lapply(fplat, rd_plat))
  if (!is.null(keep_table)) {
    info <- keep_table
    info$col <- info$gsm
    if (!is.null(pmap)) {
      hit <- info$col %in% names(pmap)
      info$col[hit] <- unname(pmap[info$col[hit]])
    }
    info <- info[info$col %in% colnames(e), ]
    info <- info[!duplicated(info$col), ]
  } else {
    info <- rd_meta(fmat[1])
    if (!is.null(keep_field)) info <- info[info[[keep_field]] %in% keep_value, ]
    info$col <- info$gsm
  }
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

raw <- read_excel("GSE8650_JIAHC.xlsx", col_names = FALSE)
sel <- data.frame(gsm = as.character(raw[[1]]),
                  illness = as.character(raw[[12]]),
                  symptoms = as.character(raw[[13]]),
                  treatment = as.character(raw[[14]]),
                  platform = as.character(raw[[24]]),
                  stringsAsFactors = FALSE)
sel <- sel[!is.na(sel$gsm) & grepl("^GSM", sel$gsm), ]
sel$illness <- sub("^[^:]+:\\s*", "", sel$illness)
sel$symptoms <- sub("^[^:]+:\\s*", "", sel$symptoms)
sel$treatment <- sub("^[^:]+:\\s*", "", sel$treatment)

g2 <- prep_gse(c("GSE8650-GPL96_series_matrix.txt.gz", "GSE8650-GPL97_series_matrix.txt.gz"),
               c("GPL96-57554.txt", "GPL97-17394.txt"),
               transform = "log2", qn = TRUE, cut = 8,
               keep_table = sel)

expr <- g2$expr
info2 <- g2$info

info2$group <- factor(ifelse(info2$illness == "Healthy", "HC", "sJIA"), c("HC", "sJIA"))
info2$drug <- ifelse(info2$group == "HC", NA, ifelse(info2$treatment == "None", "untreated", "anakinra"))

dim(expr)
table(info2$group, info2$drug)
table(info2$drug, info2$symptoms)

######2 Cell infiltration
####(1) run CIBERSORT
lin <- 2^expr - 1
write.table(round(lin, 3), "GSE8650_mix.txt", sep = "\t", quote = FALSE, col.names = NA)
source("CIBERSORT.R")
res <- CIBERSORT("LM22.txt", "GSE8650_mix.txt", perm = 100, QN = FALSE)

####(2) pairwise tests per cell type
cel_test <- function(f, grp, lev, lab) {
  d <- data.frame(cell = rep(colnames(f), each = nrow(f)),
                  val = as.vector(as.matrix(f)),
                  grp = rep(as.character(grp), times = ncol(f)),
                  stringsAsFactors = FALSE)
  d$cell <- factor(d$cell, levels = colnames(f))
  d$grp <- factor(d$grp, levels = lev, labels = lab)
  cb <- list(c(1, 2), c(1, 3), c(2, 3))
  out <- lapply(levels(d$cell), function(k) {
    s <- d[d$cell == k, ]
    do.call(rbind, lapply(seq_along(cb), function(j) {
      a <- cb[[j]][1]; b <- cb[[j]][2]
      s2 <- s[as.integer(s$grp) %in% c(a, b), ]
      s2$grp <- droplevels(s2$grp)
      p <- if (max(s$val) == 0 || nlevels(s2$grp) < 2) NA else tryCatch(wilcox.test(val ~ grp, data = s2)$p.value, error = function(e) NA)
      data.frame(cell = k, pair = paste(lab[a], "vs", lab[b]),
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
                    pt_size = 1.1, box_w = 0.6, fill_alpha = 1) {
  dk <- function(x, k = 0.62) {
    m <- col2rgb(x) / 255
    rgb(m[1] * k, m[2] * k, m[3] * k)
  }
  nn <- table(factor(as.character(grp), levels = lev))
  xl <- setNames(paste0(lab, "\n(n=", as.integer(nn[lev]), ")"), lab)
  d <- data.frame(cell = rep(colnames(f), each = nrow(f)),
                  val = as.vector(as.matrix(f)),
                  grp = rep(as.character(grp), times = ncol(f)),
                  stringsAsFactors = FALSE)
  d$cell <- factor(d$cell, levels = colnames(f))
  d$grp <- factor(d$grp, levels = lev, labels = lab)
  dc <- sapply(cols, dk)
  out <- lapply(cells, function(k) {
    s <- d[d$cell == k, ]
    t2 <- tab[tab$cell == k & !is.na(tab$star) & tab$star != "", ]
    g <- ggplot(s, aes(grp, val)) +
      geom_boxplot(aes(fill = grp), width = box_w, outlier.shape = NA, linewidth = 0.3, colour = "black") +
      geom_point(aes(colour = grp), shape = 19, size = pt_size, alpha = 0.85,
                 position = position_jitterdodge(jitter.width = box_w * 0.2, dodge.width = box_w)) +
      scale_fill_manual(values = ggplot2::alpha(cols, fill_alpha)) +
      scale_colour_manual(values = dc) +
      scale_x_discrete(labels = xl) +
      guides(colour = "none") +
      labs(x = NULL, y = "Estimated proportion", title = k) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 10, colour = "black", hjust = 0.5),
            legend.title = element_blank(),
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
  names(out) <- cells
  out
}

f <- res[, 1:22, drop = FALSE]
info2$grp3 <- factor(ifelse(info2$group == "HC", "HC", info2$drug), c("HC", "untreated", "anakinra"))

tab <- cel_test(f, info2$grp3, c("HC", "untreated", "anakinra"), c("HC", "UT", "ANA"))
tab[tab$p < 0.05, ]
tab[tab$padj < 0.05, ]

pl <- pl_cell(f, info2$grp3, c(HC = "grey90", UT = "#BB0021", ANA = "#6699CC"),
              c("HC", "untreated", "anakinra"), c("HC", "UT", "ANA"),
              c("T cells CD4 naive", "T cells CD4 memory resting"), tab,
              pt_size = 2.0, box_w = 0.42, fill_alpha = 0.65)
pl[[1]]
pl[[2]]


####3 Pathway heatmap
# load("E:/通路基因集合/gene4pathway.Rdata")
# IFN28 <- c("DDX60", "EPSTI1", "HERC5", "HERC6", "IFI27", "IFI44", "IFI44L", "IFI6",
#            "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "LAMP3", "LY6E", "MX1",
#            "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "RTP4", "SIGLEC1", "SPATS2L",
#            "USP18", "CXCL10", "GBP1", "SOCS1")
# length(intersect(IFN28, rownames(expr)))
# length(intersect(NFKB$gene, rownames(expr)))
# 
# hm <- function(e, gs, grp, cols, ttl) {
#   g <- intersect(gs, rownames(e))
#   z <- t(scale(t(e[g, , drop = FALSE])))
#   z[z > 2] <- 2
#   z[z < -2] <- -2
#   o <- order(as.integer(factor(grp, levels = names(cols))))
#   d <- data.frame(gene = rep(rownames(z), times = ncol(z)),
#                   sample = rep(colnames(z), each = nrow(z)),
#                   val = as.vector(z),
#                   grp = rep(as.character(grp), each = nrow(z)),
#                   stringsAsFactors = FALSE)
#   d$gene <- factor(d$gene, levels = rev(g))
#   d$sample <- factor(d$sample, levels = colnames(z)[o])
#   d$grp <- factor(d$grp, levels = names(cols))
#   ggplot(d, aes(sample, gene, fill = val)) +
#     geom_tile() +
#     facet_grid(~grp, scales = "free_x", space = "free_x") +
#     scale_fill_gradient2(low = "#4466B0", mid = "white", high = "#DF5734",
#                          midpoint = 0, limits = c(-2, 2), name = "z-score") +
#     labs(x = NULL, y = NULL, title = paste0(ttl, " (n=", length(g), ")")) +
#     theme_classic() +
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.line = element_blank(),
#           axis.text.y = element_text(size = 6, colour = "black"),
#           plot.title = element_text(size = 10, colour = "black", hjust = 0.5),
#           strip.background = element_blank(),
#           strip.text = element_text(size = 9, colour = "black"),
#           legend.title = element_text(size = 8),
#           legend.text = element_text(size = 7))
# }
# 
# h1 <- hm(expr, IFN28, info2$grp3, c(HC = "grey90", UT = "#DF5734", ANA = "#6699CC"), "IFN response")
# h2 <- hm(expr, NFKB$gene, info2$grp3, c(HC = "grey90", UT = "#DF5734", ANA = "#6699CC"), "NF-kB target")
# h1
# h2

####4 Differential expression
library(limma)

grp <- factor(as.character(info2$grp3), c("HC", "untreated", "anakinra"))
design <- model.matrix(~0 + grp)
colnames(design) <- levels(grp)

cm <- makeContrasts(ANA_UT = anakinra - untreated,
                    ANA_HC = anakinra - HC,
                    UT_HC  = untreated - HC,
                    levels = design)

fit <- lmFit(expr, design)
fit2 <- eBayes(contrasts.fit(fit, cm))

de <- lapply(colnames(cm), function(k) {
  tt <- topTable(fit2, coef = k, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt
})
names(de) <- colnames(cm)

####(1) summary and gene set overlap
smry <- do.call(rbind, lapply(names(de), function(k) {
  d <- de[[k]]
  s <- d[d$adj.P.Val < 0.05 & abs(d$logFC) > 1, ]
  data.frame(contrast = k,
             n_up = sum(s$logFC > 1), n_dn = sum(s$logFC < -1),
             n_p05 = sum(d$P.Value < 0.05),
             IFN28 = sum(s$gene %in% IFN28),
             IFN28_up = sum(s$gene %in% IFN28 & s$logFC > 0),
             NFKB = sum(s$gene %in% NFKB$gene),
             NFKB_up = sum(s$gene %in% NFKB$gene & s$logFC > 0),
             stringsAsFactors = FALSE)
}))
smry


#####GSE8650, sJIA PBMC, HC / untreated / anakinra
# 1 preprocessing    GPL96 + GPL97 merged, log2, quantile, low-expression filter
# 2 cell infiltration CIBERSORT, unpaired Wilcoxon, three pairwise comparisons
# 3 pathway heatmap  IFN28 and NFKB z-score
# 4 differential expr limma, reference-level contrasts vs HC
# 5 GSEA             HALLMARK + KEGG, ranked by moderated t
# 6 GSEA plots       gseaplot2, two pathways and combined

########1. preprocessing
library(limma)
library(readxl)

prep_gse <- function(fmat, fplat, transform = c("none", "log2"), qn = TRUE, cut = 5, frac = 0.2,
                     keep_field = NULL, keep_value = NULL, keep_table = NULL) {
  if (length(fmat) != length(fplat)) stop("fmat and fplat length mismatch")
  transform <- match.arg(transform)
  
  rd_expr <- function(f) {
    ln <- readLines(f)
    tb <- ln[(grep("^!series_matrix_table_begin", ln) + 1):(grep("^!series_matrix_table_end", ln) - 1)]
    con <- textConnection(tb)
    e <- as.matrix(read.delim(con, row.names = 1, check.names = FALSE))
    close(con)
    if (transform == "log2") e <- log2(e + 1)
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
  
  rd_meta <- function(f) {
    ln <- readLines(f)
    sm <- strsplit(sub("^!Sample_", "", ln[grepl("^!Sample_", ln)]), "\t", fixed = TRUE)
    k <- vapply(sm, `[`, "", 1)
    gsm <- gsub('"', "", sm[[which(k == "geo_accession")[1]]][-1])
    ch <- t(do.call(rbind, lapply(sm[which(k == "characteristics_ch1")], function(x) gsub('"', "", x[-1], fixed = TRUE))))
    colnames(ch) <- sub(":\\s.*$", "", ch[1, ])
    data.frame(gsm = gsm, ch[-1, , drop = FALSE], stringsAsFactors = FALSE)
  }
  
  e <- lapply(fmat, rd_expr)
  pmap <- NULL
  if (length(e) == 2) {
    ia <- as.integer(sub("GSM", "", colnames(e[[1]])))
    ib <- as.integer(sub("GSM", "", colnames(e[[2]])))
    j <- match(ib - 1, ia)
    if (any(is.na(j))) stop("GSM pairing failed")
    pmap <- setNames(colnames(e[[1]])[j], colnames(e[[2]]))
    colnames(e[[2]]) <- unname(pmap)
    e <- rbind(e[[1]], e[[2]])
  } else {
    e <- e[[1]]
  }
  
  m <- do.call(rbind, lapply(fplat, rd_plat))
  if (!is.null(keep_table)) {
    info <- keep_table
    info$col <- info$gsm
    if (!is.null(pmap)) {
      hit <- info$col %in% names(pmap)
      info$col[hit] <- unname(pmap[info$col[hit]])
    }
    info <- info[info$col %in% colnames(e), ]
    info <- info[!duplicated(info$col), ]
  } else {
    info <- rd_meta(fmat[1])
    if (!is.null(keep_field)) info <- info[info[[keep_field]] %in% keep_value, ]
    info$col <- info$gsm
  }
  
  e <- e[, info$col, drop = FALSE]
  
  sel <- !is.na(m$SYMBOL) & m$SYMBOL != "" & m$SYMBOL != "---" & m$PROBEID %in% rownames(e)
  m <- m[sel, ]
  m <- m[order(-rowMeans(e)[m$PROBEID]), ]
  m <- m[!duplicated(m$SYMBOL), ]
  
  e2 <- e[m$PROBEID, , drop = FALSE]
  rownames(e2) <- m$SYMBOL
  e2 <- e2[complete.cases(e2), , drop = FALSE]
  e2 <- e2[rowSums(e2 > cut) >= ceiling(frac * ncol(e2)), , drop = FALSE]
  
  list(expr = e2, info = info)
}

raw <- read_excel("GSE8650_JIAHC.xlsx", col_names = FALSE)
sel <- data.frame(gsm = as.character(raw[[1]]),
                  illness = as.character(raw[[12]]),
                  symptoms = as.character(raw[[13]]),
                  treatment = as.character(raw[[14]]),
                  platform = as.character(raw[[24]]),
                  stringsAsFactors = FALSE)
sel <- sel[!is.na(sel$gsm) & grepl("^GSM", sel$gsm), ]
sel$illness <- sub("^[^:]+:\\s*", "", sel$illness)
sel$symptoms <- sub("^[^:]+:\\s*", "", sel$symptoms)
sel$treatment <- sub("^[^:]+:\\s*", "", sel$treatment)

g2 <- prep_gse(c("GSE8650-GPL96_series_matrix.txt.gz", "GSE8650-GPL97_series_matrix.txt.gz"),
               c("GPL96-57554.txt", "GPL97-17394.txt"),
               transform = "log2", qn = TRUE, cut = 8,
               keep_table = sel)

expr <- g2$expr
info2 <- g2$info

info2$group <- factor(ifelse(info2$illness == "Healthy", "HC", "sJIA"), c("HC", "sJIA"))
info2$drug <- ifelse(info2$group == "HC", NA, ifelse(info2$treatment == "None", "untreated", "anakinra"))
info2$grp3 <- factor(ifelse(info2$group == "HC", "HC", info2$drug), c("HC", "untreated", "anakinra"))

dim(expr)
sum(is.na(expr))
table(info2$grp3)


########2. cell infiltration
####(1) deconvolution
lin <- 2^expr - 1
write.table(round(lin, 3), "GSE8650_mix.txt", sep = "\t", quote = FALSE, col.names = NA)
source("CIBERSORT.R")
res <- CIBERSORT("LM22.txt", "GSE8650_mix.txt", perm = 100, QN = FALSE)

f <- res[, 1:22, drop = FALSE]

####(2) pairwise tests per cell type
cel_test <- function(f, grp, lev, lab) {
  d <- data.frame(cell = rep(colnames(f), each = nrow(f)),
                  val = as.vector(as.matrix(f)),
                  grp = rep(as.character(grp), times = ncol(f)),
                  stringsAsFactors = FALSE)
  d$cell <- factor(d$cell, levels = colnames(f))
  d$grp <- factor(d$grp, levels = lev, labels = lab)
  cb <- list(c(1, 2), c(1, 3), c(2, 3))
  out <- lapply(levels(d$cell), function(k) {
    s <- d[d$cell == k, ]
    do.call(rbind, lapply(seq_along(cb), function(j) {
      a <- cb[[j]][1]; b <- cb[[j]][2]
      s2 <- s[as.integer(s$grp) %in% c(a, b), ]
      s2$grp <- droplevels(s2$grp)
      p <- if (max(s$val) == 0 || nlevels(s2$grp) < 2) NA else tryCatch(wilcox.test(val ~ grp, data = s2)$p.value, error = function(e) NA)
      data.frame(cell = k, pair = paste(lab[a], "vs", lab[b]),
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
                    pt_size = 1.1, box_w = 0.6, fill_alpha = 1) {
  dk <- function(x, k = 0.62) {
    m <- col2rgb(x) / 255
    rgb(m[1] * k, m[2] * k, m[3] * k)
  }
  nn <- table(factor(as.character(grp), levels = lev))
  xl <- setNames(paste0(lab, "\n(n=", as.integer(nn[lev]), ")"), lab)
  d <- data.frame(cell = rep(colnames(f), each = nrow(f)),
                  val = as.vector(as.matrix(f)),
                  grp = rep(as.character(grp), times = ncol(f)),
                  stringsAsFactors = FALSE)
  d$cell <- factor(d$cell, levels = colnames(f))
  d$grp <- factor(d$grp, levels = lev, labels = lab)
  dc <- sapply(cols, dk)
  out <- lapply(cells, function(k) {
    s <- d[d$cell == k, ]
    t2 <- tab[tab$cell == k & !is.na(tab$star) & tab$star != "", ]
    g <- ggplot(s, aes(grp, val)) +
      geom_boxplot(aes(fill = grp), width = box_w, outlier.shape = NA, linewidth = 0.3, colour = "black") +
      geom_point(aes(colour = grp), shape = 19, size = pt_size, alpha = 0.85,
                 position = position_jitterdodge(jitter.width = box_w * 0.2, dodge.width = box_w)) +
      scale_fill_manual(values = ggplot2::alpha(cols, fill_alpha)) +
      scale_colour_manual(values = dc) +
      scale_x_discrete(labels = xl) +
      guides(colour = "none") +
      labs(x = NULL, y = "Estimated proportion", title = k) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 10, colour = "black", hjust = 0.5),
            legend.title = element_blank(),
            legend.position = "top")
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
  names(out) <- cells
  out
}

####(4) run
tab <- cel_test(f, info2$grp3, c("HC", "untreated", "anakinra"), c("HC", "UT", "ANA"))
tab[tab$p < 0.05, ]
tab[tab$padj < 0.05, ]

pl <- pl_cell(f, info2$grp3, c(HC = "grey90", UT = "#DF5734", ANA = "#6699CC"),
              c("HC", "untreated", "anakinra"), c("HC", "UT", "ANA"),
              c("T cells CD4 naive", "T cells CD4 memory resting"), tab,
              pt_size = 1.0, box_w = 0.62, fill_alpha = 0.85)
pl[[1]]
pl[[2]]


########3. pathway heatmap
# load("E:/通路基因集合/gene4pathway.Rdata")
# IFN28 <- c("DDX60", "EPSTI1", "HERC5", "HERC6", "IFI27", "IFI44", "IFI44L", "IFI6",
#            "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "LAMP3", "LY6E", "MX1",
#            "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "RTP4", "SIGLEC1", "SPATS2L",
#            "USP18", "CXCL10", "GBP1", "SOCS1")
# 
# length(intersect(IFN28, rownames(expr)))
# length(intersect(NFKB$gene, rownames(expr)))
# 
# hm <- function(e, gs, grp, ttl, sz = 6) {
#   g <- intersect(gs, rownames(e))
#   z <- t(scale(t(e[g, , drop = FALSE])))
#   z[z > 2] <- 2
#   z[z < -2] <- -2
#   o <- order(as.integer(factor(grp, levels = c("HC", "untreated", "anakinra"))))
#   d <- data.frame(gene = rep(rownames(z), times = ncol(z)),
#                   sample = rep(colnames(z), each = nrow(z)),
#                   val = as.vector(z),
#                   grp = rep(as.character(grp), each = nrow(z)),
#                   stringsAsFactors = FALSE)
#   d$gene <- factor(d$gene, levels = rev(g))
#   d$sample <- factor(d$sample, levels = colnames(z)[o])
#   d$grp <- factor(d$grp, levels = c("HC", "untreated", "anakinra"))
#   ggplot(d, aes(sample, gene, fill = val)) +
#     geom_tile() +
#     facet_grid(~grp, scales = "free_x", space = "free_x") +
#     scale_fill_gradient2(low = "#4466B0", mid = "white", high = "#DF5734",
#                          midpoint = 0, limits = c(-2, 2), name = "z-score") +
#     labs(x = NULL, y = NULL, title = paste0(ttl, " (n=", length(g), ")")) +
#     theme_classic() +
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.line = element_blank(),
#           axis.text.y = element_text(size = sz, colour = "black"),
#           plot.title = element_text(size = 10, colour = "black", hjust = 0.5),
#           strip.background = element_blank(),
#           strip.text = element_text(size = 9, colour = "black"))
# }
# 
# p_ifn <- hm(expr, IFN28, info2$grp3, "IFN response")
# p_nfkb <- hm(expr, NFKB$gene, info2$grp3, "NF-kB target", sz = 3.5)
# p_ifn
# p_nfkb


########4. differential expression
grp <- factor(as.character(info2$grp3), c("HC", "untreated", "anakinra"))
design <- model.matrix(~ grp, info2)
fit <- eBayes(lmFit(expr, design))

info2$grpR <- relevel(grp, "untreated")
designR <- model.matrix(~ grpR, info2)
fitR <- eBayes(lmFit(expr, designR))

de <- list()
de$UT_HC  <- topTable(fit,  coef = "grpuntreated", number = Inf, sort.by = "P")
de$ANA_HC <- topTable(fit,  coef = "grpanakinra",  number = Inf, sort.by = "P")
de$ANA_UT <- topTable(fitR, coef = "grpRanakinra", number = Inf, sort.by = "P")

####(1) direction check
chk <- c("ISG15", "IFI44L", "MX1", "S100A8", "S100A9")
gm <- function(g, lev) sapply(g, function(x) mean(expr[x, grp == lev]))

data.frame(gene = chk,
           t_ANA_HC = de$ANA_HC$t[match(chk, rownames(de$ANA_HC))],
           d_ANA_HC = gm(chk, "anakinra") - gm(chk, "HC"),
           t_UT_HC  = de$UT_HC$t[match(chk, rownames(de$UT_HC))],
           d_UT_HC  = gm(chk, "untreated") - gm(chk, "HC"))

####(2) summary
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


########5. GSEA
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

lab <- c(ANA_UT = "higher_after_anakinra",
         ANA_HC = "higher_in_anakinra",
         UT_HC  = "higher_in_untreated")

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
table(gdp$contrast, gdp$dir)
write.xlsx(gdp,file="GSE8650_GSEAresults.xlsx")


########6. GSEA plots
library(clusterProfiler)
library(enrichplot)
gv <- function(id, nm) {
  r <- gdp[gdp$ID == id & gdp$contrast == nm, ]
  if (!nrow(r)) return(NULL)
  r[which.min(r$p.adjust), ]
}

unique(gdp$contrast)
c1 <- "UT_HC"
g1 <- gv("TNFA_SIGNALING_VIA_NFKB", c1)
g2 <- gv("INFLAMMATORY_RESPONSE", c1)
g2<-gv("INTERFERON_ALPHA_RESPONSE",c1)

p1 <- if (is.null(g1)) NULL else gseaplot2(gse[[c1]], geneSetID = "TNFA_SIGNALING_VIA_NFKB",
                                           title = sprintf("TNFA signalling via NF-kB   NES = %.2f, FDR = %.2g", g1$NES, g1$p.adjust),
                                           color = "#4466B0", pvalue_table = TRUE, base_size = 9)

p2 <- if (is.null(g2)) NULL else gseaplot2(gse[[c1]], geneSetID = "INFLAMMATORY_RESPONSE",
                                           title = sprintf("Inflammatory response   NES = %.2f, FDR = %.2g", g2$NES, g2$p.adjust),
                                           color = "#4466B0", pvalue_table = TRUE, base_size = 9)

p3 <- if (is.null(g1) || is.null(g2)) NULL else gseaplot2(gse[[c1]],
                                                          geneSetID = c("TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE"),
                                                          title = sprintf("TNFA/NF-kB (NES = %.2f, FDR = %.2g)  vs  Inflammatory (NES = %.2f, FDR = %.2g)",
                                                                          g1$NES, g1$p.adjust, g2$NES, g2$p.adjust),
                                                          color = c("#4466B0", "#FEC260"), pvalue_table = FALSE, base_size = 9)

p1
p2
p3
