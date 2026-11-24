library(Seurat)
library(SeuratObject)
library(nichenetr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(tibble)

options(future.globals.maxSize = 8 * 1024^3)

setwd("E:\\Cohort PPT\\JIA\\code\\CellCom")
pbmc1<-readRDS("E:\\Cohort PPT\\JIA\\code\\pbmc1.RDS")

##### HELPER ACCESSORS
.pick <- function(x, nm) if (is.null(x)) get(nm, envir = .GlobalEnv) else x
.getdata <- function(o) tryCatch(SeuratObject::LayerData(o, assay = "RNA", layer = "data"),
                                 error = function(e) GetAssayData(o, assay = "RNA", slot = "data"))

#####1 NICHENET RUN: receiver vs sender, patient vs control
nn_run <- function(receiver, mut, sender = NULL, n_sender = 6, obj = NULL,
                   lfc = 0.5, pct = 0.25, p_adj = 0.05,
                   top_n_ligands = 20,
                   ltm_mat = NULL, lr_net = NULL, w_net = NULL,
                   drop_pat = "^(RPL|RPS|MRPL|MRPS|MT-|HSP|DNAJ|HB[ABP]|MALAT1|NEAT1|XIST)",
                   save = FALSE, tag = NULL) {
  obj     <- if (is.null(obj)) get("pbmc1", envir = .GlobalEnv) else obj
  ltm_mat <- .pick(ltm_mat, "ltm"); lr_net <- .pick(lr_net, "lrn"); w_net <- .pick(w_net, "wn")
  mut_long <- mut
  if (is.null(sender)) {
    cc <- chat_list[[mut]]
    w  <- cc@net$weight
    s  <- sort(w[, receiver], decreasing = TRUE)
    sender <- head(names(s)[names(s) != receiver], n_sender)
  }
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
  
  out <- list(nn = nn, deg = deg, lig_expr = lig_expr, rec_expr = rec_expr,
              receiver = receiver, comparison = mut, sender = sender,
              ltm = ltm_mat, lr = lr_net,
              params = list(lfc = lfc, pct = pct, p_adj = p_adj, top_n_ligands = top_n_ligands))
  if (save) saveRDS(out, paste0("nn_", gsub("[^A-Za-z0-9]+", "_", receiver), "_", mut, ".RDS"))
  out
}

NN_DIR <- "E:/组学数据分析/细胞通讯分析学习/NicheNet"
ltm <- readRDS(file.path(NN_DIR, "ligand_target_matrix.rds"))
lrn <- readRDS(file.path(NN_DIR, "lr_network.rds"))
wn  <- readRDS(file.path(NN_DIR, "weighted_networks.rds"))
ls(pattern = "^(ltm|lrn|wn)$")

short_lab <- c("HC","Poly","S366Ter","I67R","M513K","C495Y","E250K"); names(short_lab) <- mut_levels
mut_levels <- c("Healthy controls","Polygenic JIA","S366Ter/N19KfsTer31",
                "I67R/K781TfsTer7","M513K","C495Y","E250K")

pbmc1$group_chat <- factor(pbmc1$mutation, levels = mut_levels)
table(pbmc1$group_chat, useNA = "ifany")

####(1) LPIN2, receiver cells: CD14 Monocyte; sender cells: cytotoxic CD8T, NK 
x1 <- nn_run(
  receiver = "CD14 Monocyte",
  mut      = "I67R/K781TfsTer7",
  sender   = c("Cytotoxic CD8 T"),
  top_n_ligands = 60)
x2 <- nn_run(
  receiver = "CD14 Monocyte",
  mut      = "I67R/K781TfsTer7",
  sender   = c("NK"),
  top_n_ligands = 60)

saveRDS(x1,file="LPIN2_CD14_CD8T_Niche.RDS")
saveRDS(x2,file="LPIN2_CD14_NK_Niche.RDS")

####(2) PSTPIP1, 
#sender cells: CD16 Monocyte; receiver cells: CD14 Monocyte
#sender cells: pDC; receiver cells: CD14 Monocyte 
y1 <- nn_run(
  receiver = "CD14 Monocyte",
  mut      = "E250K",
  sender   = c("CD16 Monocyte"),
  top_n_ligands = 60)
y2 <- nn_run(
  receiver = "CD14 Monocyte",
  mut      = "E250K",
  sender   = c("pDC"),
  top_n_ligands = 60)

saveRDS(y1,file="PSTPIP1_CD14_CD16_Niche.RDS")
saveRDS(y2,file="PSTPIP1_CD14_pDC_Niche.RDS")

#####2 plot function for nichenet analysis
nn_fig <- function(x, ligands = NULL, genes = NULL, top_n_lig = 20, top_n_gene = 60,
                   row_cut = 1e-4, col_cut = 1e-4,
                   w = c(3.0, 7.5, 11), h = c(4.2, 4.4), tile_w = 1, tile_h = 1,
                   col_lig = "#6699CC", col_target = "#2B4B8F", col_lr = "#AF478A",
                   col_low = "whitesmoke", x_from = NULL, base_size = 8,
                   title = NULL, file = NULL, width = 15, height = 9) {
  la <- x$nn$ligand_activities
  la <- la[order(-la$pearson), ]
  k  <- min(top_n_lig, nrow(la))
  la <- la[seq_len(k), ]
  la$test_ligand <- factor(la$test_ligand, levels = rev(la$test_ligand))
  pA <- ggplot(la, aes(y = test_ligand)) +
    geom_segment(aes(x = if (is.null(x_from)) 0 else x_from, xend = pearson, yend = test_ligand),
                 colour = "grey85", linewidth = 0.45) +
    geom_point(aes(x = pearson), size = 2.1, colour = col_lig) +
    scale_x_continuous(limits = c(x_from, NA), expand = expansion(mult = c(0.01, 0.06))) +
    theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"), axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)) +
    labs(x = "Ligand activity (Pearson)", y = NULL)
  
  ligB <- if (is.null(ligands)) as.character(la$test_ligand)[seq_len(min(10, k))] else ligands
  ligB <- intersect(ligB, rownames(x$ltm))
  deg_g <- intersect(x$deg$gene, colnames(x$ltm))
  if (is.null(genes)) {
    sc  <- apply(x$ltm[ligB, deg_g, drop = FALSE], 2, max, na.rm = TRUE)
    genes <- names(sort(sc, decreasing = TRUE))[seq_len(min(top_n_gene, length(sc)))]
  } else genes <- intersect(genes, deg_g)
  dfB <- expand.grid(ligand = ligB, target = genes, stringsAsFactors = FALSE)
  dfB$rp <- x$ltm[cbind(dfB$ligand, dfB$target)]
  
  ##### filter B by row/col sum thresholds
  matB <- matrix(NA_real_, nrow = length(ligB), ncol = length(genes),
                 dimnames = list(ligB, genes))
  matB[cbind(dfB$ligand, dfB$target)] <- dfB$rp
  keep_lB <- names(which(rowSums(matB, na.rm = TRUE) >= row_cut))
  keep_gB <- names(which(colSums(matB, na.rm = TRUE) >= col_cut))
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
  
  ##### filter C by row/col sum thresholds
  matC <- xtabs(pot ~ ligand + receptor, data = lrm)
  keep_lC <- names(which(rowSums(matC) >= row_cut))
  keep_rC <- names(which(colSums(matC) >= col_cut))
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


objs <- list(x1=x1, x2=x2, y1=y1, y2=y2)
tags <- c(x1="LPIN2_CD14_CD8T", x2="LPIN2_CD14_NK",
          y1="PSTPIP1_CD14_CD16", y2="PSTPIP1_CD14_pDC")
titles <- c(x1="LPIN2 | CD14 Monocyte <- Cytotoxic CD8 T", x2="LPIN2 | CD14 Monocyte <- NK",
            y1="PSTPIP1 | CD14 Monocyte <- CD16 Monocyte", y2="PSTPIP1 | CD14 Monocyte <- pDC")

fig <- list()
for (nm in names(objs)) {
  fig[[nm]] <- nn_fig(objs[[nm]], row_cut=0.01, col_cut=0.005, title=titles[nm])
  print(fig[[nm]])
  d <- attr(fig[[nm]], "data")
  openxlsx::write.xlsx(list(LigandActivity=d$A, LigandTarget=d$B, LigandReceptor=d$C),
                       file=paste0(tags[nm], "_Niche_data.xlsx"))
}
fig$x1 | fig$x2; fig$y1 | fig$y2






# 
# ##################################################################################
# ##################################################################################
# #####total  pipline
# library(Seurat); library(nichenetr); library(tidyverse)
# library(clusterProfiler)
# setwd("E:/Cohort PPT/JIA/code/CellCom")
# 
# ct_keep <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
#              "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC")
# mut_levels <- c("Healthy controls","Polygenic JIA","S366Ter/N19KfsTer31",
#                 "I67R/K781TfsTer7","M513K","C495Y","E250K")
# short_lab <- c("HC","Poly","S366Ter","I67R","M513K","C495Y","E250K")
# names(short_lab) <- mut_levels
# sel <- c("MHC-I","MHC-II","ICAM","ITGB2","SELPLG","CD99","LCK","CCL","CLEC","GALECTIN",
#          "IFN-II","ADGRE5","TNF","IL1")
# 
# NN_DIR <- "E:/组学数据分析/细胞通讯分析学习/NicheNet"
# ltm <- readRDS(file.path(NN_DIR, "ligand_target_matrix.rds"))
# lrn <- readRDS(file.path(NN_DIR, "lr_network.rds"))
# wn  <- readRDS(file.path(NN_DIR, "weighted_networks.rds"))
# 
# pbmc1 <- readRDS("E:/Cohort PPT/JIA/code/pbmc1.RDS")
# pbmc1 <- subset(pbmc1, celltype %in% ct_keep)
# pbmc1$celltype <- factor(pbmc1$celltype, levels = ct_keep)
# if (!"group_chat" %in% colnames(pbmc1@meta.data)) {
#   pbmc1$group_chat <- factor(pbmc1$mutation, levels = mut_levels)
# }
# pbmc1$group_chat <- factor(pbmc1$group_chat, levels = mut_levels)
# DefaultAssay(pbmc1) <- "RNA"
# 
# chat_list <- lapply(short_lab, function(g) readRDS(paste0("chat_", g, ".RDS")))
# names(chat_list) <- short_lab
# 
# incoming <- sapply(chat_list, function(cc) colSums(cc@net$weight))
# outgoing <- sapply(chat_list, function(cc) rowSums(cc@net$weight))
# rank_receiver <- sort(rowSums(incoming), decreasing = TRUE)
# rank_sender   <- sort(rowSums(outgoing), decreasing = TRUE)
# print(round(rank_receiver, 3)); print(round(rank_sender, 3))
# 
# sender_from_cc <- function(mut, receiver, n = 6) {
#   cc <- chat_list[[mut]]
#   if (is.null(cc)) cc <- chat_list$HC
#   w <- cc@net$weight
#   if (!receiver %in% colnames(w)) return(setdiff(ct_keep, receiver))
#   s <- sort(w[, receiver], decreasing = TRUE)
#   s <- s[names(s) != receiver]
#   head(names(s), n)
# }
# 
# run_nn <- function(receiver, mut, tag, sender, top_n = 20, lfc = 0.15, pct = 0.10) {
#   mut_long <- names(short_lab)[short_lab == mut]
#   sub <- subset(pbmc1, subset = group_chat %in% c("Healthy controls", mut_long) &
#                   celltype %in% c(receiver, sender))
#   sub$status <- factor(ifelse(sub$group_chat == "Healthy controls", "control", "patient"),
#                        levels = c("control", "patient"))
#   Idents(sub) <- "celltype"
#   
#   print(table(sub$celltype, sub$status))
#   print(round(prop.table(table(sub$celltype, sub$status), 2), 3))
#   if (min(table(sub$status)) < 10) {
#     message("  skipped: <10 cells in one condition")
#     return(NULL)
#   }
#   
#   nn <- nichenet_seuratobj_aggregate(
#     seurat_obj = sub, top_n_ligands = top_n, receiver = receiver, sender = sender,
#     condition_colname = "status", condition_oi = "patient", condition_reference = "control",
#     ligand_target_matrix = ltm, lr_network = lrn, weighted_networks = wn,
#     assay_oi = "RNA", expression_pct = pct, lfc_cutoff = lfc)
#   
#   nn$.receiver <- receiver
#   nn$.comparison <- mut
#   saveRDS(nn, paste0("nn_", tag, ".RDS"))
#   
#   p1 <- nn$ligand_activities %>% dplyr::top_n(top_n, pearson) %>%
#     arrange(desc(pearson)) %>% mutate(test_ligand = factor(test_ligand, levels = test_ligand)) %>%
#     ggplot(aes(test_ligand, pearson)) + geom_col(fill = "#6699CC") +
#     coord_flip() + labs(x = "Ligand", y = "Pearson", title = tag) +
#     theme_minimal(base_size = 9)
#   
#   pdf(paste0("nn_", tag, ".pdf"), width = 12, height = 9)
#   print(p1)
#   if (!is.null(nn$ligand_target_heatmap))
#     print(nn$ligand_target_heatmap +
#             scale_fill_gradient2(low = "whitesmoke", high = "royalblue", name = "Regulatory potential"))
#   if (!is.null(nn$ligand_receptor_heatmap)) print(nn$ligand_receptor_heatmap)
#   dev.off()
#   
#   nn
# }
# 
# receivers   <- names(head(rank_receiver, 8))
# comparisons <- short_lab[short_lab != "HC"]
# scan <- list()
# 
# for (mut in comparisons) {
#   for (rv in receivers) {
#     tag <- paste0(gsub("[^A-Za-z0-9]+", "_", rv), "_", mut)
#     message("[", format(Sys.time(), "%H:%M:%S"), "] ", tag)
#     sender <- sender_from_cc(mut, rv, n = 6)
#     scan[[tag]] <- tryCatch(run_nn(rv, mut, tag, sender),
#                             error = function(e) list(error = conditionMessage(e)))
#     gc()
#   }
# }
# 
# summ <- map_dfr(names(scan), function(k) {
#   x <- scan[[k]]
#   if (is.null(x$ligand_activities)) return(tibble())
#   x$ligand_activities %>% dplyr::top_n(10, pearson) %>%
#     mutate(receiver = x$.receiver, comparison = x$.comparison)
# })
# print(dim(summ))
# 
# inter  <- CellChatDB.human$interaction
# lig2pw <- inter %>% dplyr::select(ligand, pathway_name) %>%
#   filter(pathway_name %in% sel) %>% distinct(ligand, pathway_name)
# 
# summ2 <- summ %>% left_join(lig2pw, by = c("test_ligand" = "ligand"))
# summ2 <- summ2[order(summ2$comparison, summ2$receiver, -summ2$pearson), ]
# print(head(summ2, 20))
# write.csv(summ2, "nn_scan_top_ligands.csv", row.names = FALSE)
# saveRDS(scan, "nn_scan_all.RDS")
