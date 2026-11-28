library(Seurat); library(SeuratObject); library(nichenetr)
library(ggplot2); library(patchwork); library(Matrix); library(dplyr)
library(tibble); library(openxlsx)

options(future.globals.maxSize = 8 * 1024^3)
setwd("E:\\Cohort PPT\\JIA\\code\\CellCom")
pbmc1 <- readRDS("E:\\Cohort PPT\\JIA\\code\\pbmc1.RDS")

mut_levels <- c("Healthy controls", "Polygenic JIA", "S366Ter/N19KfsTer31",
                "I67R/K781TfsTer7", "M513K", "C495Y", "E250K")
short_lab <- c("HC", "Poly", "S366Ter", "I67R", "M513K", "C495Y", "E250K")
names(short_lab) <- mut_levels
pbmc1$group_chat <- factor(pbmc1$mutation, levels = mut_levels)

NN_DIR <- "E:/组学数据分析/细胞通讯分析学习/NicheNet/V2"
ltm <- readRDS(file.path(NN_DIR, "ligand_target_matrix_nsga2r_final.rds"))
lrn <- readRDS(file.path(NN_DIR, "lr_network_human_21122021.rds"))
wn  <- readRDS(file.path(NN_DIR, "weighted_networks_nsga2r_final.rds"))

ifn_lig <- "IFNG"
ifn_rec <- c("IFNGR1", "IFNGR2")
ifn_tg  <- c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2", "OAS1",
             "STAT1", "IRF1", "IRF7", "CXCL10")
il1_lig <- c("IL1A", "IL1B", "IL1RN")
il1_rec <- c("IL1R1", "IL1R2", "IL1RAP")
il1_tg  <- c("IL6", "CXCL8", "S100A8", "S100A9")

#####1 nn_run: receiver vs sender, patient vs control
.pick <- function(x, nm) if (is.null(x)) get(nm, envir = .GlobalEnv) else x
.getdata <- function(o) tryCatch(SeuratObject::LayerData(o, assay = "RNA", layer = "data"),
                                 error = function(e) GetAssayData(o, assay = "RNA", slot = "data"))

nn_run <- function(receiver, mut, sender = NULL, n_sender = 6, obj = NULL,
                   lfc = 0.1, pct = 0.05, p_adj = 0.05,
                   top_n_ligands = 60,
                   ltm_mat = NULL, lr_net = NULL, w_net = NULL,
                   drop_pat = "^(RPL|RPS|MRPL|MRPS|MT-|HSP|DNAJ|HB[ABP]|MALAT1|NEAT1|XIST)") {
  obj     <- if (is.null(obj)) get("pbmc1", envir = .GlobalEnv) else obj
  ltm_mat <- .pick(ltm_mat, "ltm"); lr_net <- .pick(lr_net, "lrn"); w_net <- .pick(w_net, "wn")
  mut_long <- mut
  
  sub <- subset(obj, subset = group_chat %in% c("Healthy controls", mut_long) &
                  celltype %in% c(receiver, sender))
  sub$status <- factor(ifelse(sub$group_chat == "Healthy controls", "control", "patient"),
                       levels = c("control", "patient"))
  if (min(table(sub$status)) < 10) stop("fewer than 10 cells in one condition")
  
  mat <- .getdata(sub)
  snd_cells <- colnames(sub)[sub$celltype %in% sender]
  rec_cells <- colnames(sub)[sub$celltype == receiver]
  lig_expr <- Matrix::rowMeans(mat[, snd_cells, drop = FALSE])
  rec_expr <- Matrix::rowMeans(mat[, rec_cells, drop = FALSE])
  
  Idents(sub) <- "celltype"
  nn <- nichenet_seuratobj_aggregate(
    seurat_obj = sub, top_n_ligands = top_n_ligands, receiver = receiver, sender = sender,
    condition_colname = "status", condition_oi = "patient", condition_reference = "control",
    ligand_target_matrix = ltm_mat, lr_network = lr_net, weighted_networks = w_net,
    assay_oi = "RNA", expression_pct = pct, lfc_cutoff = lfc)
  
  sub_r <- sub[, rec_cells]
  Idents(sub_r) <- "status"
  deg <- FindMarkers(sub_r, ident.1 = "patient", ident.2 = "control",
                     logfc.threshold = 0, verbose = FALSE)
  deg$gene <- rownames(deg)
  deg <- deg[deg$p_val_adj < p_adj & abs(deg$avg_log2FC) >= lfc &
               pmax(deg$pct.1, deg$pct.2) >= pct, ]
  deg <- deg[!grepl(drop_pat, deg$gene), ]
  deg <- deg[order(-abs(deg$avg_log2FC)), ]
  
  list(nn = nn, deg = deg, lig_expr = lig_expr, rec_expr = rec_expr,
       receiver = receiver, comparison = mut, sender = sender,
       ltm = ltm_mat, lr = lr_net,
       params = list(lfc = lfc, pct = pct, p_adj = p_adj, top_n_ligands = top_n_ligands))
}

##### 2 Runing Nichenet analysis for LPIN2 AND PSTPIP1
###(1) LPIN2: CD14 Monocyte <- Cytotoxic CD8 T / NK, pathway IFN-II
x1 <- nn_run(receiver = "CD14 Monocyte", mut = "I67R/K781TfsTer7",
             sender = c("Cytotoxic CD8 T"), lfc = 0.05, pct = 0.25, p_adj = 0.05,
             top_n_ligands = 60)
saveRDS(x1, file = "LPIN2_CD14_CD8T_Niche.RDS")

x2 <- nn_run(receiver = "CD14 Monocyte", mut = "I67R/K781TfsTer7",
             sender = c("NK"), lfc = 0.05, pct = 0.25, p_adj = 0.05,
             top_n_ligands = 60)
saveRDS(x2, file = "LPIN2_CD14_NK_Niche.RDS")

###(2)PSTPIP1: CD14 Monocyte <- CD16 Monocyte / pDC, pathway IL1
y1 <- nn_run(receiver = "CD14 Monocyte", mut = "E250K",
             sender = c("CD16 Monocyte"), lfc = 0.05, pct = 0.25, p_adj = 0.05,
             top_n_ligands = 60)
saveRDS(y1, file = "PSTPIP1_CD14_CD16_Niche.RDS")

y2 <- nn_run(receiver = "CD14 Monocyte", mut = "E250K",
             sender = c("pDC"), lfc = 0.05, pct = 0.25, p_adj = 0.05,
             top_n_ligands = 60)
saveRDS(y2, file = "PSTPIP1_CD14_pDC_Niche.RDS")

#####3 plot function, retain=FALSE threshold only; retain=TRUE force keep pathway genes
nn_fig <- function(x, ligands = NULL, genes = NULL, top_n_lig = 20, top_n_gene = 60,
                   lt_row_cut = 1e-4, lt_col_cut = 1e-4,
                   lr_row_cut = 1e-4, lr_col_cut = 1e-4,
                   keep_lig = NULL, keep_rec = NULL, keep_tg = NULL,
                   retain = FALSE,
                   w = c(3.0, 7.5, 11), h = c(4.2, 4.4), tile_w = 1, tile_h = 1,
                   col_lig = "#007ABA", col_target = "#AF478A", col_lr = "#AF478A",
                   col_keep = "#EC706E", col_low = "whitesmoke", x_from = NULL,
                   base_size = 8, title = NULL, file = NULL, width = 15, height = 9) {
  
  if (!isTRUE(retain)) { keep_lig <- NULL; keep_rec <- NULL; keep_tg <- NULL }
  
  la_all <- x$nn$ligand_activities
  la_all <- la_all[order(-la_all$pearson), ]
  top_lig <- as.character(la_all$test_ligand)[seq_len(min(top_n_lig, nrow(la_all)))]
  la <- la_all[la_all$test_ligand %in% top_lig, ]
  la$retained <- FALSE
  if (!is.null(keep_lig)) {
    add <- setdiff(intersect(keep_lig, la_all$test_ligand), la$test_ligand)
    if (length(add)) {
      la_add <- la_all[la_all$test_ligand %in% add, ]
      la_add$retained <- TRUE
      la <- rbind(la, la_add)
    }
  }
  la <- la[order(-la$pearson), ]
  la$test_ligand <- factor(la$test_ligand, levels = rev(la$test_ligand))
  
  pA <- ggplot(la, aes(y = test_ligand)) +
    geom_segment(aes(x = if (is.null(x_from)) 0 else x_from, xend = pearson, yend = test_ligand),
                 colour = "grey85", linewidth = 0.45) +
    geom_point(aes(x = pearson, colour = retained), size = 2.1) +
    scale_colour_manual(values = c("FALSE" = col_lig, "TRUE" = col_keep), guide = "none") +
    scale_x_continuous(limits = c(x_from, NA), expand = expansion(mult = c(0.01, 0.06))) +
    theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"), axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) +
    labs(x = "Ligand activity (Pearson)", y = NULL)
  
  ligB <- if (is.null(ligands)) as.character(la$test_ligand)[seq_len(min(10, nrow(la)))] else ligands
  ligB <- intersect(ligB, colnames(x$ltm))
  if (!is.null(keep_lig)) ligB <- union(ligB, intersect(keep_lig, colnames(x$ltm)))
  deg_g <- intersect(x$deg$gene, rownames(x$ltm))
  if (is.null(genes)) {
    sc  <- apply(x$ltm[deg_g, ligB, drop = FALSE], 1, max, na.rm = TRUE)
    genes <- names(sort(sc, decreasing = TRUE))[seq_len(min(top_n_gene, length(sc)))]
  } else genes <- intersect(genes, deg_g)
  if (!is.null(keep_tg)) genes <- union(genes, intersect(keep_tg, deg_g))
  
  dfB <- expand.grid(ligand = ligB, target = genes, stringsAsFactors = FALSE)
  dfB$rp <- x$ltm[cbind(dfB$target, dfB$ligand)]
  dfB$retained <- dfB$ligand %in% keep_lig | dfB$target %in% keep_tg
  
  matB <- matrix(NA_real_, nrow = length(ligB), ncol = length(genes),
                 dimnames = list(ligB, genes))
  matB[cbind(dfB$ligand, dfB$target)] <- dfB$rp
  keep_lB <- names(which(rowSums(matB, na.rm = TRUE) >= lt_row_cut))
  keep_gB <- names(which(colSums(matB, na.rm = TRUE) >= lt_col_cut))
  if (!is.null(keep_lig)) keep_lB <- union(keep_lB, intersect(keep_lig, ligB))
  if (!is.null(keep_tg))  keep_gB <- union(keep_gB, intersect(keep_tg, genes))
  dfB <- dfB[dfB$ligand %in% keep_lB & dfB$target %in% keep_gB, ]
  
  pB <- ggplot(dfB, aes(x = target, y = ligand, fill = rp)) +
    geom_tile(width = tile_w, height = tile_h) +
    scale_fill_gradient(low = col_low, high = col_target, name = "Regulatory\npotential") +
    scale_y_discrete(limits = rev(keep_lB)) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size - 1),
          axis.text.y = element_text(colour = "black"), panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) +
    labs(x = "Predicted target genes", y = NULL)
  
  lr <- x$lr[as.character(x$lr$from) %in% ligB, c("from", "to")]
  lr <- lr[!duplicated(paste(lr$from, lr$to)), ]
  rmean <- function(r) {
    v <- sapply(strsplit(as.character(r), "_"), function(z) {
      e <- x$rec_expr[z]; if (any(is.na(e))) NA_real_ else min(e) })
    as.numeric(v)
  }
  lrm <- data.frame(ligand = as.character(lr$from), receptor = as.character(lr$to),
                    lig_expr = as.numeric(x$lig_expr[as.character(lr$from)]),
                    rec_expr = rmean(lr$to), stringsAsFactors = FALSE)
  lrm <- lrm[is.finite(lrm$lig_expr) & is.finite(lrm$rec_expr), ]
  lrm$pot <- (lrm$lig_expr / max(lrm$lig_expr)) * (lrm$rec_expr / max(lrm$rec_expr))
  lrm <- lrm[lrm$pot > 0, ]
  lrm$retained <- lrm$ligand %in% keep_lig | lrm$receptor %in% keep_rec
  
  matC <- xtabs(pot ~ ligand + receptor, data = lrm)
  keep_lC <- names(which(rowSums(matC) >= lr_row_cut))
  keep_rC <- names(which(colSums(matC) >= lr_col_cut))
  if (!is.null(keep_lig)) keep_lC <- union(keep_lC, intersect(keep_lig, unique(lrm$ligand)))
  if (!is.null(keep_rec)) keep_rC <- union(keep_rC, intersect(keep_rec, unique(lrm$receptor)))
  lrm <- lrm[lrm$ligand %in% keep_lC & lrm$receptor %in% keep_rC, ]
  
  pC <- ggplot(lrm, aes(x = receptor, y = ligand, fill = pot)) +
    geom_tile(width = tile_w, height = tile_h) +
    scale_fill_gradient(low = col_low, high = col_lr, name = "Interaction\npotential",
                        limits = c(0, 1)) +
    scale_y_discrete(limits = rev(intersect(keep_lC, unique(lrm$ligand)))) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size - 1),
          axis.text.y = element_text(colour = "black"), panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) +
    labs(x = "Receptors", y = NULL)
  
  p <- (pA | pB) / pC + plot_layout(widths = w, heights = h)
  if (!is.null(title))
    p <- p + plot_annotation(title = title, theme = theme(plot.title = element_text(size = base_size + 2)))
  
  attr(p, "parts") <- list(A = pA, B = pB, C = pC)
  attr(p, "data")  <- list(A = la, B = dfB, C = lrm)
  if (!is.null(file)) ggsave(file, p, width = width, height = height)
  p
}

#####4 two modes: raw threshold only, and retain forced pathway genes
objs <- list(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
tags <- c(x1 = "LPIN2_CD14_CD8T", x2 = "LPIN2_CD14_NK",
          y1 = "PSTPIP1_CD14_CD16", y2 = "PSTPIP1_CD14_pDC")
titles <- c(x1 = "LPIN2 | CD14 Monocyte <- Cytotoxic CD8 T",
            x2 = "LPIN2 | CD14 Monocyte <- NK",
            y1 = "PSTPIP1 | CD14 Monocyte <- CD16 Monocyte",
            y2 = "PSTPIP1 | CD14 Monocyte <- pDC")
keep_map <- list(
  x1 = list(lig = ifn_lig, rec = ifn_rec, tg = ifn_tg),
  x2 = list(lig = ifn_lig, rec = ifn_rec, tg = ifn_tg),
  y1 = list(lig = il1_lig, rec = il1_rec, tg = il1_tg),
  y2 = list(lig = il1_lig, rec = il1_rec, tg = il1_tg))

fig_raw <- list()
for (nm in names(objs)) {
  km <- keep_map[[nm]]
  fig_raw[[nm]] <- nn_fig(objs[[nm]], retain = FALSE,
                          lt_row_cut = 1e-4, lt_col_cut = 1e-4,
                          lr_row_cut = 1e-4, lr_col_cut = 1e-4,
                          keep_lig = km$lig, keep_rec = km$rec, keep_tg = km$tg,
                          title = paste(titles[nm], "| raw"))
  print(fig_raw[[nm]])
  d <- attr(fig_raw[[nm]], "data")
  openxlsx::write.xlsx(list(LigandActivity = d$A, LigandTarget = d$B, LigandReceptor = d$C),
                       file = paste0(tags[nm], "_raw_Niche_data.xlsx"))
}
fig_raw$x1 | fig_raw$x2; fig_raw$y1 | fig_raw$y2

fig_ret <- list()
for (nm in names(objs)) {
  km <- keep_map[[nm]]
  fig_ret[[nm]] <- nn_fig(objs[[nm]], retain = TRUE,
                          lt_row_cut = 1e-4, lt_col_cut = 1e-4,
                          lr_row_cut = 1e-4, lr_col_cut = 1e-4,
                          keep_lig = km$lig, keep_rec = km$rec, keep_tg = km$tg,
                          title = paste(titles[nm], "| retain"))
  print(fig_ret[[nm]])
  d <- attr(fig_ret[[nm]], "data")
  openxlsx::write.xlsx(list(LigandActivity = d$A, LigandTarget = d$B, LigandReceptor = d$C),
                       file = paste0(tags[nm], "_retain_Niche_data.xlsx"))
}
fig_ret$x1 | fig_ret$x2
fig_ret$y1 | fig_ret$y2

######5 return each single figure (whole figure and its A/B/C panels)
fig_all <- list(raw = fig_raw, retain = fig_ret)
for (md in names(fig_all)) {
  for (nm in names(fig_all[[md]])) {
    p <- fig_all[[md]][[nm]]
    parts <- attr(p, "parts")
    cat("\n---", md, tags[nm], "---\n")
    print(p)          # whole figure
    print(parts$A)    # ligand activity
    print(parts$B)    # ligand-target regulatory potential
    print(parts$C)    # ligand-receptor interaction potential
  }
}

get_part <- function(p, which = c("A","B","C")) attr(p, "parts")[[match.arg(which)]]
p1<-get_part(fig_ret$x1, "B")
p2<-get_part(fig_ret$x2, "B")

p3<-get_part(fig_ret$y1, "B")
p4<-get_part(fig_ret$y2, "B")

(p1 | p2) / (p3 | p4) + plot_layout(guides = "collect")

#######6 check gene expressions, produce violin plot for batch
library(Seurat)
pbmc1<-readRDS("E:/Cohort PPT/JIA/code/pbmc1.RDS")
IFN_gene<-c("ADAR","B4GALT5","BMPR2","BTN3A1","BTN3A2","CASP1","CCL5","CD47","CEBPD","CHST11","CITED2","DDX58","DUSP6","ERAP2","FOS","GADD45B","GALNT10","GBP1","GBP2","HLA−A","HLA−C","HLA−DRA","HLA−E","ICAM1","IFI16","IFITM3","IL15","IRF1","IRF2","ISG15","JAK2","JUN","JUNB","KLF4","KLF6","KLF7","LAP3","MAP3K5","MICAL2","MX1","NEDD9","PFKFB3","PLAUR","PLSCR1","PMAIP1","PSMB8","PSME1","PSTPIP2","RBCK1","SECTM1","SMURF2","SOD2","SP100","SP110","STAT1","STAT3","STK38L","TBC1D8","TGFBI","TNFAIP2","TNFSF10","VCAN")
IL_gene<-c("ADAM19","AHR","ALOX5AP","AMPD3","BHLHE40","BID","BIRC2","BMPR2","CD82","CFLAR","CHST11","DDIT4","ELK3","FOS","HMOX1","ICAM1","IER3","IFNGR2","IRF1","ITGA5","ITGAV","IVNS1ABP","KLF10","KLF7","LIMK2","LYN","MAP3K5","MICAL2","NAMPT","NEDD9","NFKB1","NFKB2","NFKBIA","NFKBIZ","NINJ1","PALLD","PDLIM7","PFKFB3","PLAUR","S100A8","S100A9","SAT1","SKIL","SLC25A37","SMAD3","SMAD7","SMURF2","STAT5A","STK38L","SUSD6","TBC1D8","TGFBI","TGFBR1","TNFAIP2","TNFAIP3","TNFAIP8","TNIP1","TNS1","ULK1","VDR","VEGFA","WTAP")

#####(1) gene sets and contrast mapping
library(Seurat); library(ggplot2); library(patchwork)

IFN_gene <- c("ADAR","B4GALT5","BMPR2","BTN3A1","BTN3A2","CASP1","CCL5","CD47","CEBPD",
              "CHST11","CITED2","DDX58","DUSP6","ERAP2","FOS","GADD45B","GALNT10","GBP1","GBP2",
              "HLA-A","HLA-C","HLA-DRA","HLA-E","ICAM1","IFI16","IFITM3","IL15","IRF1","IRF2",
              "ISG15","JAK2","JUN","JUNB","KLF4","KLF6","KLF7","LAP3","MAP3K5","MICAL2","MX1",
              "NEDD9","PFKFB3","PLAUR","PLSCR1","PMAIP1","PSMB8","PSME1","PSTPIP2","RBCK1",
              "SECTM1","SMURF2","SOD2","SP100","SP110","STAT1","STAT3","STK38L","TBC1D8","TGFBI",
              "TNFAIP2","TNFSF10","VCAN")

IL1_gene <- c("ADAM19","AHR","ALOX5AP","AMPD3","BHLHE40","BID","BIRC2","BMPR2","CD82",
              "CFLAR","CHST11","DDIT4","ELK3","FOS","HMOX1","ICAM1","IER3","IFNGR2","IRF1",
              "ITGA5","ITGAV","IVNS1ABP","KLF10","KLF7","LIMK2","LYN","MAP3K5","MICAL2","NAMPT",
              "NEDD9","NFKB1","NFKB2","NFKBIA","NFKBIZ","NINJ1","PALLD","PDLIM7","PFKFB3","PLAUR",
              "S100A8","S100A9","SAT1","SKIL","SLC25A37","SMAD3","SMAD7","SMURF2","STAT5A",
              "STK38L","SUSD6","TBC1D8","TGFBI","TGFBR1","TNFAIP2","TNFAIP3","TNFAIP8","TNIP1",
              "TNS1","ULK1","VDR","VEGFA","WTAP")

IFN_gene <- gsub("−", "-", IFN_gene)
IL1_gene <- gsub("−", "-", IL1_gene)

mut_hc  <- "Healthy controls"
mut_i67 <- "I67R/K781TfsTer7"
mut_e25 <- "E250K"
cell_use <- "CD14 Monocyte"

#####(2) Seurat VlnPlot function: HC vs one mutation, CD14 Monocyte
vio_seurat <- function(sce, genes, mut_pat, pat_lab,
                       col_hc = "#6699CC", col_pat = "#E2A9C9",
                       cell = cell_use, assay = "RNA") {
  genes <- genes[genes %in% rownames(sce[[assay]])]
  if (length(genes) == 0) stop("no valid genes in assay")
  
  md <- sce@meta.data
  keep <- rownames(md)[md$celltype == cell & md$mutation %in% c(mut_hc, mut_pat)]
  sub <- sce[, keep]
  md2 <- sub@meta.data
  md2$plot_mut <- factor(ifelse(as.character(md2$mutation) == mut_hc, "HC", pat_lab),
                         levels = c("HC", pat_lab))
  sub@meta.data <- md2
  
  n_hc  <- sum(md2$plot_mut == "HC")
  n_pat <- sum(md2$plot_mut == pat_lab)
  cat("== contrast", pat_lab, "| HC n:", n_hc, "|", pat_lab, "n:", n_pat, "\n")
  if (n_pat == 0) warning(pat_lab, " cells not matched, check mutation column")
  
  fill_vals <- setNames(c(col_hc, col_pat), c("HC", pat_lab))
  expr <- tryCatch(GetAssayData(sub, assay = assay, layer = "data"),
                   error = function(e) GetAssayData(sub, assay = assay, slot = "data"))
  plots <- list(); stats_rows <- list()
  
  for (g in genes) {
    v <- as.numeric(expr[g, ])
    ps <- as.character(md2$plot_mut)
    
    p <- VlnPlot(sub, features = g, group.by = "plot_mut", assay = assay, pt.size = 0) +
      scale_fill_manual(values = fill_vals) +
      theme_classic() + xlab(NULL) + ylab(NULL) + ggtitle(g) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.ticks.x = element_blank(),
            plot.margin = margin(2, 4, 2, 2))
    plots[[g]] <- p
    
    for (gp in c("HC", pat_lab)) {
      vv <- v[ps == gp]
      if (length(vv) < 2) next
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        gene = g, celltype = cell, contrast = pat_lab, group = gp,
        n = length(vv), mean_expr = mean(vv), median_expr = median(vv),
        pct_expressing = mean(vv > 0) * 100, stringsAsFactors = FALSE)
    }
  }
  list(plots = plots, stats = do.call(rbind, stats_rows))
}

#####(3) save pdf: 8 rows x 10 columns per page, enlarged panels
save_vio_pdf <- function(plots, pdf_file, ncol = 10, nrow = 8,
                         col_width = 1.6, row_height = 1.6) {
  if (length(plots) == 0) { warning("no plots to save"); return(invisible(NULL)) }
  dir.create(dirname(pdf_file), showWarnings = FALSE, recursive = TRUE)
  per_page <- ncol * nrow
  grp_idx <- split(seq_along(plots), ceiling(seq_along(plots) / per_page))
  pdf(pdf_file, width = ncol * col_width, height = nrow * row_height)
  on.exit(dev.off(), add = TRUE)
  for (i in grp_idx) print(wrap_plots(plots[i], ncol = ncol))
  cat("pdf saved:", pdf_file, "| genes:", length(plots),
      "| pages:", length(grp_idx), "| exists:", file.exists(pdf_file), "\n")
}

out_dir <- "E:/Cohort PPT/JIA/code/CellPreprocess"
res_ifn <- vio_seurat(pbmc1, IFN_gene, mut_i67, "I67R", col_pat = "#E2A9C9")
res_il1 <- vio_seurat(pbmc1, IL1_gene, mut_e25, "E250K", col_pat = "#7B1FA2")

save_vio_pdf(res_ifn$plots, file.path(out_dir, "IFN_gene_CD14_HC_vs_I67R.pdf"))
save_vio_pdf(res_il1$plots, file.path(out_dir, "IL1_gene_CD14_HC_vs_E250K.pdf"))

print(res_ifn$stats, row.names = FALSE)
print(res_il1$stats, row.names = FALSE)

#####(4) showing gene sets and mapping
IFN_show <- c("DDX58","GBP1","JAK2","MX1","FOS","JUNB","KLF4",
              "ERAP2","GADD45B","KLF6","SMURF2")

IL1_show <- c("ALOX5AP","ICAM1","ITGAV","ITGA5","VEGFA","BHLHE40","HMOX1",
              "BIRC2","CD82","IER3","SMAD3","STAT5A","ULK1","VDR",
              "KLF10","SMAD7","SMURF2","DDIT4","ELK3","BMPR2")

mut_hc  <- "Healthy controls"
mut_i67 <- "I67R/K781TfsTer7"
mut_e25 <- "E250K"
cell_use <- "CD14 Monocyte"

library(Seurat); library(ggplot2); library(patchwork)

vio_panel <- function(sce, genes, mut_pat, pat_lab, group_colors, ncol,
                      cell = cell_use, assay = "RNA") {
  genes <- genes[genes %in% rownames(sce[[assay]])]
  if (length(genes) == 0) stop("no valid genes in assay")
  
  md <- sce@meta.data
  keep <- rownames(md)[md$celltype == cell & md$mutation %in% c(mut_hc, mut_pat)]
  sub <- sce[, keep]
  md2 <- sub@meta.data
  md2$plot_mut <- factor(ifelse(as.character(md2$mutation) == mut_hc, "HC", pat_lab),
                         levels = c("HC", pat_lab))
  sub@meta.data <- md2
  
  n_hc  <- sum(md2$plot_mut == "HC")
  n_pat <- sum(md2$plot_mut == pat_lab)
  cat("== contrast", pat_lab, "| HC n:", n_hc, "|", pat_lab, "n:", n_pat, "\n")
  if (n_pat == 0) warning(pat_lab, " cells not matched, check mutation column")
  
  fill_vals <- group_colors[c("HC", pat_lab)]
  expr <- tryCatch(GetAssayData(sub, assay = assay, layer = "data"),
                   error = function(e) GetAssayData(sub, assay = assay, slot = "data"))
  plots <- list(); stats_rows <- list()
  
  for (g in genes) {
    v <- as.numeric(expr[g, ])
    ps <- as.character(md2$plot_mut)
    
    p <- VlnPlot(sub, features = g, group.by = "plot_mut", assay = assay, pt.size = 0) +
      scale_fill_manual(values = fill_vals) +
      theme_classic() + xlab(NULL) + ylab(NULL) + ggtitle(g) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 10, colour = "black"),
            plot.margin = margin(2, 4, 2, 2))
    plots[[g]] <- p
    
    for (gp in c("HC", pat_lab)) {
      vv <- v[ps == gp]
      if (length(vv) < 2) next
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        gene = g, celltype = cell, contrast = pat_lab, group = gp,
        n = length(vv), mean_expr = mean(vv), median_expr = median(vv),
        pct_expressing = mean(vv > 0) * 100, stringsAsFactors = FALSE)
    }
  }
  panel <- wrap_plots(plots, ncol = ncol)
  list(panel = panel, plots = plots, stats = do.call(rbind, stats_rows))
}

col_ifn <- c("HC" = "#6699CC", "I67R" = "#E2A9C9")
col_il1 <- c("HC" = "#6699CC", "E250K" = "#FFEFC1")

res_ifn <- vio_panel(pbmc1, IFN_show, mut_i67, "I67R", col_ifn, ncol = length(IFN_show))
res_il1 <- vio_panel(pbmc1, IL1_show, mut_e25, "E250K", col_il1, ncol = 10)

res_ifn$panel
res_il1$panel

print(res_ifn$stats, row.names = FALSE)
print(res_il1$stats, row.names = FALSE)

