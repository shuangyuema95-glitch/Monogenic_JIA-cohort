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
                   col_lig = "#6699CC", col_target = "#6699CC", col_lr = "#AF478A",
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
p1<-get_part(fig_raw$x1, "B")
p2<-get_part(fig_raw$x2, "B")

p3<-get_part(fig_raw$y1, "B")
p4<-get_part(fig_raw$y1, "B")
