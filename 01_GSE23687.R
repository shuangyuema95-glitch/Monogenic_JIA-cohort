#####GSE23687, JIA patients with MTX treatment
# 1 preprocessing        series matrix -> symbol matrix, RMA log2, 20% filter
# 2 cell infiltration    CIBERSORT on 2^x linear, paired + group Wilcoxon
# 3 pathway heatmap      IFN28 and NFKB z-score, 4 subject columns
# 4 differential expr    limma, post-pre paired, R70-NR at month 6
# 5 GSEA                 HALLMARK + KEGG, ranked by moderated t, FDR < 0.05
# 6 GSEA plots           gseaplot2, per pathway and combined
# contrast direction     post_vs_pre = timepost (post - pre, 11 pairs)
#                        R_vs_NR_post = groupR70 (R70 - NR, 7 vs 4)


dir <- "E:/Cohort PPT/JIA/bulkdata"; setwd(dir)

########1. preprocessing
####(1) read expression matrix
ln <- readLines("GSE23687_series_matrix.txt.gz")
tb <- ln[(grep("^!series_matrix_table_begin", ln) + 1):(grep("^!series_matrix_table_end", ln) - 1)]
con <- textConnection(tb); e1 <- as.matrix(read.delim(con, row.names = 1, check.names = FALSE)); close(con)

####(2) sample annotation
smp <- strsplit(sub("^!Sample_", "", ln[grepl("^!Sample_", ln)]), "\t", fixed = TRUE)
k <- vapply(smp, `[`, "", 1)
gsm <- gsub('"', "", smp[[which(k == "geo_accession")[1]]][-1])
ch <- t(do.call(rbind, lapply(smp[which(k == "characteristics_ch1")], function(x) gsub('"', "", x[-1], fixed = TRUE))))
colnames(ch) <- sub(":\\s.*$", "", ch[1, ]); rownames(ch) <- gsm
info <- data.frame(gsm = gsm, subject = ch[, "subject"], resp = ch[, "responder"], tx = ch[, "treatment"], stringsAsFactors = FALSE)
info$time <- factor(ifelse(grepl("pre", info$tx, ignore.case = TRUE), "pre", "post"), c("pre", "post"))
info$group <- ifelse(grepl("not", info$resp, ignore.case = TRUE), "NR", "R70")
rownames(info) <- gsm

####(3) platform annotation
pl <- readLines("GPL570-55999.txt")
pb <- grep("^!platform_table_begin", pl); pe <- grep("^!platform_table_end", pl)
pt <- if (length(pb)) pl[(pb + 1):(pe - 1)] else pl[grep("^ID\t", pl)[1]:length(pl)]
ph <- strsplit(pt[1], "\t", fixed = TRUE)[[1]]; sx <- grep("gene symbol", ph, ignore.case = TRUE)[1]

####(4) probe to symbol
x2 <- strsplit(pt[-1], "\t", fixed = TRUE)
s2 <- lapply(x2, function(z) trimws(strsplit(z[sx], "///", fixed = TRUE)[[1]]))
m <- data.frame(PROBEID = rep(vapply(x2, `[`, "", 1), lengths(s2)), SYMBOL = trimws(unlist(s2)), stringsAsFactors = FALSE)
m <- m[!is.na(m$SYMBOL) & m$SYMBOL != "" & m$SYMBOL != "---" & m$PROBEID %in% rownames(e1), ]
m <- m[order(-rowMeans(e1)[m$PROBEID]), ]; m <- m[!duplicated(m$SYMBOL), ]

####(5) symbol-level matrix
e2 <- e1[m$PROBEID, ]; rownames(e2) <- m$SYMBOL
e2f <- e2[rowSums(e2 > 5) >= ceiling(0.2 * ncol(e2)), ]

####(6) align and save
ord <- order(info$subject, info$time)
e2 <- e2[, ord]; e2f <- e2f[, ord]; info <- info[ord, ]
saveRDS(list(expr = e2f, expr_all = e2, info = info), "GSE23687_pre.rds")


#######2. cell infiltration
####(1) deconvolution
# e2 is RMA log2; CIBERSORT needs non-log linear input
write.table(round(2^e2, 3), "GSE23687_mix.txt", sep = "\t", quote = FALSE, col.names = NA)
source("CIBERSORT.R")
frac <- CIBERSORT("LM22.txt", "GSE23687_mix.txt", perm = 0, QN = TRUE)
frac <- frac[info$gsm, ]

####(2) two comparisons
g <- function(s, t) info$gsm[match(paste(s, t), paste(info$subject, info$time))]
s70 <- unique(info$subject[info$group == "R70"])
s40 <- unique(info$subject[info$group == "NR"])
pre <- g(c(s70, s40), "pre"); post <- g(c(s70, s40), "post")
q70 <- g(s70, "post"); q40 <- g(s40, "post")
j <- 1:22
pt <- vapply(j, function(k) suppressWarnings(wilcox.test(frac[post, k], frac[pre, k], paired = TRUE)$p.value), 0)
pg <- vapply(j, function(k) suppressWarnings(wilcox.test(frac[q70, k], frac[q40, k], exact = FALSE)$p.value), 0)
tab <- data.frame(cell = colnames(frac)[j],
                  m_pre = colMeans(frac[pre, j]), m_post = colMeans(frac[post, j]),
                  p_time = pt, padj_time = p.adjust(pt, "BH"),
                  m_R70 = colMeans(frac[q70, j]), m_NR = colMeans(frac[q40, j]),
                  p_grp = pg, padj_grp = p.adjust(pg, "BH"), check.names = FALSE)
write.csv(tab, "GSE23687_cell_test.csv", row.names = FALSE)

####(3) plot function
pl2 <- function(d, cols, pv, paired = FALSE, ylab = "Cell proportion (%)") {
  t1 <- tapply(d$pct, d$cell, max)
  r1 <- tapply(d$pct, d$cell, function(v) diff(range(v)))
  r1[r1 == 0] <- 1
  lab <- data.frame(cell = names(t1), p = as.numeric(pv[names(t1)]), stringsAsFactors = FALSE)
  lab$star <- as.character(cut(lab$p, c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns")))
  lab$star[is.na(lab$p)] <- "NA"
  lab$y <- as.numeric(t1[lab$cell]) + as.numeric(r1[lab$cell]) * 0.18
  lab$yend <- lab$y - as.numeric(r1[lab$cell]) * 0.05
  p <- ggplot(d, aes(xx, pct, fill = ff)) +
    geom_boxplot(width = 0.55, outlier.shape = NA)
  if (paired) p <- p + geom_line(aes(group = subject), colour = "grey60", linewidth = 0.3)
  p + geom_point(size = 0.9, position = position_jitter(width = 0.06)) +
    geom_segment(data = lab, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
    geom_segment(data = lab, aes(x = 1, xend = 1, y = yend, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
    geom_segment(data = lab, aes(x = 2, xend = 2, y = yend, yend = y), inherit.aes = FALSE, linewidth = 0.35) +
    geom_text(data = lab, aes(x = 1.5, y = y, label = star), inherit.aes = FALSE, size = 2.8, vjust = -0.3) +
    facet_wrap(~cell, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.30))) +
    coord_cartesian(ylim = c(0, NA), clip = "off") +
    labs(x = NULL, y = ylab) + theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black"), strip.background = element_blank())
}

####(4) two figures
lg <- data.frame(cell = rep(colnames(frac)[j], each = nrow(frac)), gsm = rep(rownames(frac), 22),
                 pct = as.vector(as.matrix(frac[, j])) * 100)
lg$group <- info$group[match(lg$gsm, info$gsm)]
lg$time <- info$time[match(lg$gsm, info$gsm)]
lg$subject <- info$subject[match(lg$gsm, info$gsm)]

lg1 <- lg[lg$time == "post", ]
lg1$xx <- factor(lg1$group, c("R70", "NR")); lg1$ff <- lg1$xx
p_grp <- pl2(lg1, c(R70 = "#6699CC", NR = "#DF5734"), setNames(tab$padj_grp, tab$cell))

lg2 <- lg
lg2$xx <- factor(lg2$time, c("pre", "post")); lg2$ff <- lg2$xx
p_time <- pl2(lg2, c(pre = "grey85", post = "#6699CC"), setNames(tab$padj_time, tab$cell), paired = TRUE)

saveRDS(list(frac = frac, tab = tab, p_grp = p_grp, p_time = p_time), "GSE23687_cibersort.rds")

#######3. pathway heatmap
load("E:/通路基因集合/gene4pathway.Rdata")
IFN28 <- c("DDX60", "EPSTI1", "HERC5", "HERC6", "IFI27", "IFI44", "IFI44L", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "LAMP3", "LY6E", "MX1", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "RTP4", "SIGLEC1", "SPATS2L", "USP18", "CXCL10", "GBP1", "SOCS1")
NFKB <- NFKB$gene
####(1) heatmap function
library(ggplot2)
hm <- function(gs, ttl, sz = 4) {
  g <- intersect(gs, rownames(e2))
  z <- t(scale(t(e2[g, ])))
  o <- g[hclust(dist(z))$order]
  hz <- data.frame(gene = rep(o, ncol(e2)), gsm = rep(colnames(e2), each = length(o)), z = as.vector(z[o, ]))
  hz$subject <- factor(info$subject[match(hz$gsm, info$gsm)], unique(info$subject))
  hz$group <- factor(info$group[match(hz$gsm, info$gsm)], c("R70", "NR"))
  hz$time <- factor(info$time[match(hz$gsm, info$gsm)], c("pre", "post"))
  hz$gene <- factor(hz$gene, rev(o))
  ggplot(hz, aes(subject, gene, fill = z)) + geom_tile() +
    facet_grid(. ~ group + time, scales = "free", space = "free") +
    scale_fill_gradient2(low = "#4466B0", mid = "white", high = "#DF5734", midpoint = 0) +
    labs(x = NULL, y = NULL, fill = "z-score", title = ttl) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 6),
          axis.text.y = element_text(colour = "black", size = sz),
          axis.line = element_blank(), axis.ticks = element_blank(), strip.background = element_blank())
}

p_ifn <- hm(IFN28, "IFN28")
p_nfkb <- hm(NFKB, "NFKB", sz = 3.5)

#######4. differential expression
####(1) paired fit within a subject set
library(limma)

d <- info
d$time <- factor(d$time, c("pre", "post"))
mm <- model.matrix(~ subject + time, d)
fit <- eBayes(lmFit(e2f[, d$gsm], mm))
de_pre_post <- topTable(fit, coef = "timepost", number = Inf, sort.by = "P")

dp <- d[d$time == "post", ]
dp$group <- factor(dp$group, c("NR", "R70"))
mg <- model.matrix(~ group, dp)
fg <- eBayes(lmFit(e2f[, dp$gsm], mg))
de_R_NR <- topTable(fg, coef = "groupR70", number = Inf, sort.by = "P")

de <- list(pre_vs_post_all11 = de_pre_post, R_vs_NR_post = de_R_NR)

#######5. KEGG over-representation
# library(clusterProfiler)
# library(org.Hs.eg.db)
# ora <- function(df, cut_p = 0.05, cut_lfc = 1) {
#   df <- df[!is.na(df$adj.P.Val) & df$adj.P.Val < cut_p & abs(df$logFC) > cut_lfc, ]
#   up <- rownames(df)[df$logFC > 0]
#   dn <- rownames(df)[df$logFC < 0]
#   eg <- function(g) {
#     if (length(g) < 5) return(NULL)
#     enrichKEGG(bitr(g, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID, organism = "hsa",
#                pvalueCutoff = 0.05, qvalueCutoff = 0.2)
#   }
#   list(n_up = length(up), n_down = length(dn), up = eg(up), down = eg(dn))
# }
# 
# kk <- lapply(de, ora)
# sapply(kk, function(x) c(n_up = x$n_up, n_down = x$n_down))
# 
# nm <- c(MTX_all = "postVSpre_all11", MTX_R70 = "postVSpre_R70", MTX_NR = "postVSpre_NR",
#         pre_70vs40 = "R70VS_NR_pre", post_70vs40 = "R70VS_NR_post")
# names(de) <- nm[names(de)]
# 
# kk <- lapply(de, ora)
# 
# pk <- do.call(rbind, lapply(names(kk), function(n) {
#   do.call(rbind, lapply(c("up", "down"), function(d) {
#     x <- kk[[n]][[d]]
#     if (is.null(x)) return(NULL)
#     df <- as.data.frame(setReadable(x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
#     data.frame(contrast = n, direction = d, ID = df$ID, Description = df$Description,
#                Count = df$Count, GeneRatio = df$GeneRatio, BgRatio = df$BgRatio,
#                pvalue = df$pvalue, p.adjust = df$p.adjust, geneID = df$geneID, row.names = NULL)
#   }))
# }))
# pk <- pk[order(pk$contrast, pk$pvalue), ]
# pk

######6 GSEA analysis
####(1) running enrichment
library(msigdbr)
library(clusterProfiler)

h <- msigdbr(species = "Homo sapiens", category = "H")
k <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
sets <- rbind(data.frame(term = sub("^HALLMARK_", "", h$gs_name), gene = h$gene_symbol),
              data.frame(term = sub("^KEGG_", "", k$gs_name), gene = k$gene_symbol))

gse <- lapply(de, function(x) {
  v <- setNames(x[["t"]], rownames(x))
  GSEA(sort(v[!is.na(v)], decreasing = TRUE), TERM2GENE = sets, pvalueCutoff = 0.2, verbose = FALSE)
})

lab <- c(pre_vs_post_all11 = "higher_after_MTX", R_vs_NR_post = "higher_in_R")
saveRDS(gse,file="E:/Cohort PPT/JIA/bulkdata/1112_gse.RDS")

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


####(2)GSEA plots
gv <- function(id) {
  r <- gdp[gdp$ID == id, ]
  r[which.max(r$p.adjust == min(r$p.adjust)), ]
}

g1 <- gv("TNFA_SIGNALING_VIA_NFKB")
g2 <- gv("INFLAMMATORY_RESPONSE")

p1 <- gseaplot2(gse$pre_vs_post_all11, geneSetID = "TNFA_SIGNALING_VIA_NFKB",
                title = sprintf("TNFA signalling via NF-kB   NES = %.2f, FDR = %.2g", g1$NES, g1$p.adjust),
                color = "#4466B0", pvalue_table = TRUE, base_size = 9)

p2 <- gseaplot2(gse$pre_vs_post_all11, geneSetID = "INFLAMMATORY_RESPONSE",
                title = sprintf("Inflammatory response   NES = %.2f, FDR = %.2g", g2$NES, g2$p.adjust),
                color = "#4466B0", pvalue_table = TRUE, base_size = 9)

p3 <- gseaplot2(gse$pre_vs_post_all11,
                geneSetID = c("TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE"),
                title = sprintf("TNFA/NF-kB (NES = %.2f, FDR = %.2g)  vs  Inflammatory (NES = %.2f, FDR = %.2g)",
                                g1$NES, g1$p.adjust, g2$NES, g2$p.adjust),
                color = c("#4466B0", "#FEC260"), pvalue_table = FALSE, base_size = 9)

p1
p2
p3
