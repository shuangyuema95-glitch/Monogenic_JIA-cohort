library(Seurat); library(data.table)
pbmc1 <- readRDS("pbmc1.RDS")

ct_keep <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
             "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC")
mut_levels <- c("Healthy controls","Polygenic JIA","S366Ter/N19KfsTer31",
                "I67R/K781TfsTer7","M513K","C495Y","E250K")
short_lab <- c("HC","Poly","S366Ter","I67R","M513K","C495Y","E250K"); names(short_lab) <- mut_levels

pbmc1 <- subset(pbmc1, celltype %in% ct_keep)
pbmc1$celltype   <- factor(pbmc1$celltype, levels = ct_keep)
pbmc1$group_chat <- factor(pbmc1$mutation, levels = mut_levels)
pbmc1 <- pbmc1[, which(!is.na(pbmc1$group_chat))]

for (g in mut_levels) {
  message("[", format(Sys.time(), "%H:%M:%S"), "] ", short_lab[g])
  cells_g <- colnames(pbmc1)[which(pbmc1$group_chat == g)]
  sub  <- pbmc1[, cells_g]
  cnt  <- GetAssayData(sub, assay = "RNA", slot = "counts")
  cnt  <- cnt[Matrix::rowSums(cnt > 0) >= 0.1 * ncol(cnt), ]
  meta <- data.frame(Cell = colnames(sub), cell_type = as.character(sub$celltype))
  fwrite(meta, paste0("cpdb_", short_lab[g], "_meta.txt"), sep = "\t", encoding = "UTF-8")
  m <- as.matrix(cnt)
  fwrite(data.frame(Gene = rownames(m), m, check.names = FALSE),
         paste0("cpdb_", short_lab[g], "_counts.txt"), sep = "\t", encoding = "UTF-8")
  message("  genes: ", nrow(m), " | cells: ", ncol(m))
  rm(sub, cnt, m); gc()
}
