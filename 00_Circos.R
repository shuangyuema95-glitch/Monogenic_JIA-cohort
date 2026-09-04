# ============================================
# JIA Figure 1: Circos (chr-grouped + expression preference + 8 cell tracks)
# Inside->outside: chr block -> Mb axis -> gene names -> preference -> 8 cell expr
# ============================================
setwd("E:\\Cohort PPT\\JIA")
# ---------- 0. CONFIG ----------
excel_file      <- "E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx"
hpa_sc_type     <- "E:\\Cohort PPT\\JIA\\rna_single_cell_type.tsv"
hpa_sc_cluster  <- "E:\\Cohort PPT\\JIA\\rna_single_cell_cluster.tsv"
output_pdf      <- "JIA_Figure1_circos2.pdf"
mech_colors <- c(
  "deficiency"          = "#6699CC",
  "gain-of-function"    = "#FEC260",
  "haploinsufficiency"  = "#AF478A",
  "na"                  = "grey60"
)
pref_colors <- c(
  "Myeloid-biased"  = "#BB0021",
  "Lymphoid-biased" = "#B5BBE3",
  "Ubiquitous"      = "gray90"
)
gene_alias <- c("STING1"="TMEM173")
cell_order <- c("NK","CD4T","CD8T","Bcell","Monocyte","Macrophage","Neutrophil","DC")
lymphoid_cells <- c("NK","CD4T","CD8T","Bcell")
myeloid_cells  <- c("Monocyte","Macrophage","Neutrophil","DC")
win_sector <- 5000000
# ---------- 1. PACKAGES ----------
for (pkg in c("openxlsx","circlize","ComplexHeatmap","data.table","dplyr","tidyr","tibble","grid")) {
  if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
    install.packages(pkg, repos="https://cloud.r-project.org")
    library(pkg, character.only=TRUE)
  }
}
for (pkg in c("org.Hs.eg.db","AnnotationDbi","GenomicFeatures","TxDb.Hsapiens.UCSC.hg38.knownGene")) {
  if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
    if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
    library(pkg, character.only=TRUE)
  }
}
select <- dplyr::select
# ---------- 2. READ EXCEL ----------
cat("Reading Excel:", excel_file, "\n")
summary_data <- read.xlsx(excel_file)
summary_data$gene      <- toupper(trimws(summary_data$gene))
summary_data$mechanism <- tolower(trimws(summary_data$mechanism))
summary_data$mechanism[is.na(summary_data$mechanism)] <- "na"
gene_meta <- summary_data %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(Mechanism = paste(unique(mechanism[mechanism != ""]), collapse="/"),
                   n_cases = n(), .groups="drop") %>%
  dplyr::mutate(Gene=as.character(gene),
                Mechanism=ifelse(Mechanism==""|is.na(Mechanism),"na",Mechanism)) %>%
  dplyr::select(Gene, Mechanism, n_cases)
gene_meta$Gene <- as.character(gene_meta$Gene)
gene_meta$Mechanism <- as.character(gene_meta$Mechanism)
gene_meta$line_color <- sapply(gene_meta$Mechanism, function(m)
  ifelse(m %in% names(mech_colors), mech_colors[m], "grey60"))
# ---------- 3. hg38 COORDINATES via TxDb ----------
query_genes <- unique(c(gene_meta$Gene, unlist(gene_alias)))
sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys=query_genes,
                                 columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
alias_to_canonical <- setNames(names(gene_alias), gene_alias)
sym_map$Canonical <- ifelse(sym_map$SYMBOL %in% names(alias_to_canonical),
                            alias_to_canonical[sym_map$SYMBOL], as.character(sym_map$SYMBOL))
sym_map <- sym_map[!is.na(sym_map$ENTREZID) & sym_map$ENTREZID!="", ]
sym_map <- sym_map[!duplicated(sym_map$Canonical), ]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
g <- suppressMessages(GenomicFeatures::genes(txdb))
coord_df <- as.data.frame(g)
coord_df$ENTREZID <- as.character(coord_df$gene_id)
coord_df$chr <- as.character(coord_df$seqnames)
valid_chr <- paste0("chr", c(as.character(1:22), "X", "Y"))
coord_df <- coord_df[coord_df$chr %in% valid_chr, c("ENTREZID","chr","start","end")]
gene_info <- merge(sym_map, coord_df, by="ENTREZID", all.x=TRUE)
gene_info <- gene_info[!is.na(gene_info$start), ]
gene_info <- gene_info[!duplicated(gene_info$Canonical), ]
gene_meta$chr <- NA_character_; gene_meta$start <- NA_integer_; gene_meta$end <- NA_integer_
for (i in seq_len(nrow(gene_meta))) {
  hit <- gene_info[gene_info$Canonical == gene_meta$Gene[i], ]
  if (nrow(hit) > 0) {
    gene_meta$chr[i] <- hit$chr[1]; gene_meta$start[i] <- as.integer(hit$start[1]); gene_meta$end[i] <- as.integer(hit$end[1])
  }
}
gene_meta <- gene_meta[!is.na(gene_meta$start), ]
chr_order <- paste0("chr", c(as.character(1:22), "X", "Y"))
gene_meta <- gene_meta[order(match(gene_meta$chr, chr_order), gene_meta$start), ]
rownames(gene_meta) <- NULL
cat("Genes sorted:", nrow(gene_meta), "\n")
# ---------- 4. HPA SINGLE-CELL EXPRESSION ----------
all_search_genes <- unique(c(gene_meta$Gene, unlist(gene_alias)))
process_hpa_sc <- function(fp, sg) {
  if (!file.exists(fp)) return(NULL)
  dt <- fread(fp, sep="\t", header=TRUE)
  if ("Gene name" %in% colnames(dt)) dt <- dt %>% dplyr::rename(Gene_name=`Gene name`)
  if ("Cell type" %in% colnames(dt)) dt <- dt %>% dplyr::rename(Immune_cell=`Cell type`)
  else if ("Immune cell" %in% colnames(dt)) dt <- dt %>% dplyr::rename(Immune_cell=`Immune cell`)
  dt %>% dplyr::filter(Gene_name %in% sg) %>% dplyr::mutate(nTPM=as.numeric(nTPM))
}
dt_type <- process_hpa_sc(hpa_sc_type, all_search_genes)
dt_cluster <- process_hpa_sc(hpa_sc_cluster, all_search_genes)
map_c <- function(dt) {
  if (is.null(dt)) return(NULL)
  dt$Gene_canonical <- as.character(dt$Gene_name)
  for (cn in names(gene_alias)) dt$Gene_canonical[dt$Gene_name == gene_alias[cn]] <- cn
  dt
}
dt_type <- map_c(dt_type); dt_cluster <- map_c(dt_cluster)
grp <- function(dt) {
  if (is.null(dt)) return(NULL)
  dt %>% dplyr::mutate(Cell_group=dplyr::case_when(
    Immune_cell %in% c("memory CD4 T-cell","naive CD4 T-cell","T-cells","CD4 T-cells") ~ "CD4T",
    Immune_cell %in% c("memory CD8 T-cell","naive CD8 T-cell","CD8 T-cells") ~ "CD8T",
    Immune_cell %in% c("memory B-cell","naive B-cell","B-cells","B-cell") ~ "Bcell",
    Immune_cell %in% c("NK-cell","NK-cells","NK cells") ~ "NK",
    Immune_cell %in% c("classical monocyte","intermediate monocyte","non-classical monocyte","monocytes","Monocytes") ~ "Monocyte",
    Immune_cell %in% c("myeloid DC","plasmacytoid DC","dendritic cells","Dendritic cells","mDC","pDC") ~ "DC",
    Immune_cell %in% c("neutrophil","Neutrophil","granulocytes","Granulocytes") ~ "Neutrophil",
    Immune_cell %in% c("Macrophages","macrophage","Macrophage","tissue macrophage") ~ "Macrophage",
    TRUE ~ NA_character_
  )) %>% dplyr::filter(!is.na(Cell_group))
}
dt_type <- grp(dt_type); dt_cluster <- grp(dt_cluster)
dt_list <- list(dt_type, dt_cluster); dt_list <- dt_list[!sapply(dt_list, is.null)]
dt_all <- dplyr::bind_rows(dt_list) %>%
  dplyr::distinct(Gene_canonical, Cell_group, Immune_cell, .keep_all=TRUE) %>%
  dplyr::group_by(Gene_canonical, Cell_group) %>%
  dplyr::summarise(nTPM_mean=mean(nTPM, na.rm=TRUE), .groups="drop")
mat_sc <- dt_all %>%
  tidyr::pivot_wider(names_from=Cell_group, values_from=nTPM_mean, values_fill=0) %>%
  tibble::column_to_rownames("Gene_canonical") %>% as.data.frame()
for (ct in cell_order) if (!ct %in% colnames(mat_sc)) mat_sc[[ct]] <- 0
for (g in gene_meta$Gene) if (!g %in% rownames(mat_sc)) mat_sc[g,] <- 0
mat_sc <- mat_sc[gene_meta$Gene, cell_order, drop=FALSE]
mat_sc[is.na(mat_sc)] <- 0
mat_sc_log <- log2(as.matrix(mat_sc) + 1)
# ---------- 4b. EXPRESSION PREFERENCE ----------
lymphoid_mean <- rowMeans(mat_sc_log[, lymphoid_cells])
myeloid_mean  <- rowMeans(mat_sc_log[, myeloid_cells])
pref_score <- (myeloid_mean - lymphoid_mean) / (myeloid_mean + lymphoid_mean)
pref_score[is.nan(pref_score)] <- 0
pref_class <- ifelse(pref_score > 0.15, "Myeloid-biased",
                     ifelse(pref_score < -0.15, "Lymphoid-biased", "Ubiquitous"))
gene_meta$pref_score <- pref_score
gene_meta$pref_class <- pref_class
gene_meta$pref_color <- pref_colors[pref_class]
cat("\n=== Expression Preference ===\n")
print(data.frame(Gene=gene_meta$Gene, Lymphoid=round(lymphoid_mean,2),
                 Myeloid=round(myeloid_mean,2), Score=round(pref_score,3), Class=pref_class))
cat("==============================\n\n")
# ---------- 5. BUILD SECTORS (one per gene, no blank gaps) ----------
gene_windows <- data.frame(
  chr   = gene_meta$chr,
  start = pmax(gene_meta$start - win_sector, 0),
  end   = gene_meta$end   + win_sector,
  stringsAsFactors = FALSE
)
# Within same chromosome: gap=0 (physically adjacent)
# Between different chromosomes: gap=6
gaps <- c()
for (i in 1:(nrow(gene_meta)-1)) {
  if (gene_meta$chr[i] == gene_meta$chr[i+1]) gaps <- c(gaps, 0)
  else gaps <- c(gaps, 6)
}
gaps <- c(gaps, 8)
chr_first <- c(); chr_last <- c()
for (cn in unique(gene_meta$chr)) {
  idx <- which(gene_meta$chr == cn)
  chr_first[cn] <- min(idx); chr_last[cn] <- max(idx)
}
max_sc <- max(mat_sc_log, na.rm=TRUE); if (max_sc <= 0) max_sc <- 1
expr_col_fun <- colorRamp2(breaks=c(0, max_sc*0.5, max_sc), colors=c("#007ABA", "#FFFFFF", "#BB0021"))

# ---------- 6. DRAW CIRCOS ----------
circos.clear()
pdf(output_pdf, width=14, height=14)
circos.par(cell.padding = c(0, 0, 0, 0),
           track.margin = c(0, 0),
           gap.after = gaps)
fa <- paste0(gene_windows$chr, ":", gene_windows$start, "-", gene_windows$end)
circos.initialize(factors = factor(fa, levels = fa),
                  xlim = cbind(gene_windows$start, gene_windows$end))

# --- OUTERMOST: 8 cell expression tracks ---
for (ct in rev(cell_order)) {
  local({
    ct_local <- ct
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      i <- CELL_META$sector.numeric.index
      v <- as.numeric(mat_sc_log[i, ct_local])
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col=expr_col_fun(v), border="black", lwd=0.5)
      if (i == 1) {
        circos.text(CELL_META$xlim[1], 0.5, ct_local, cex=0.45, font=1,
                    facing="inside", adj=c(0,0.5), col="black")
      }
    }, track.height = 0.04, bg.border = NA)
  })
}

# --- EXPRESSION PREFERENCE track ---
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  i <- CELL_META$sector.numeric.index
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
              col=gene_meta$pref_color[i], border="black", lwd=0.5)
  short <- ifelse(gene_meta$pref_class[i]=="Myeloid-biased","M",
                  ifelse(gene_meta$pref_class[i]=="Lymphoid-biased","L","U"))
  circos.text(mean(CELL_META$xlim), 0.5, short, cex=0.45, font=2,
              facing="inside", col="white")
}, track.height = 0.04, bg.border = NA)

# --- Gene names track ---
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  i <- CELL_META$sector.numeric.index
  circos.text(mean(CELL_META$xlim), 0.5, gene_meta$Gene[i],
              cex = 0.9, font = 1, facing = "inside",
              col = gene_meta$line_color[i])
}, track.height = 0.15, bg.border = NA)

# --- Mb coordinate axis track ---
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  x0 <- CELL_META$xlim[1]; x1 <- CELL_META$xlim[2]; width <- x1 - x0
  if (width <= 8e6) by <- 1e6
  else if (width <= 15e6) by <- 2e6
  else by <- 5e6
  start_pos <- ceiling(x0 / by) * by
  tick_n <- 0
  for (pos in seq(start_pos, x1, by)) {
    if (pos > x1) break
    circos.segments(pos, 0.6, pos, 1, lwd = 0.8)
    tick_n <- tick_n + 1
    if (tick_n %% 2 == 1) {
      mb_val <- pos / 1e6
      if (mb_val == round(mb_val)) label_text <- paste0(mb_val, "Mb")
      else label_text <- paste0(round(mb_val, 1), "Mb")
      circos.text(pos, 0.35, label_text, cex=0.4,
                  facing="inside", adj=c(0.5,0.5))
    }
  }
}, track.height = 0.12, bg.border = NA)

# --- INNERMOST: chromosome block (one continuous block, zero-gap + overlap to kill seams) ---
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  i <- CELL_META$sector.numeric.index
  cn <- gene_meta$chr[i]; fi <- chr_first[cn]; li <- chr_last[cn]
  
  # Extend rect by 0.5% on each side so adjacent same-chr rects OVERLAP,
  # completely eliminating PDF sub-pixel white seams between sectors.
  w <- CELL_META$xlim[2] - CELL_META$xlim[1]
  ov <- w * 0.005
  rect_x1 <- CELL_META$xlim[1] - ov
  rect_x2 <- CELL_META$xlim[2] + ov
  
  # Filled rectangle (no border)
  circos.rect(rect_x1, 0.2, rect_x2, 0.8, col = "grey85", border = NA)
  
  # Top and bottom borders on every sector (continuous across zero-gap)
  circos.segments(rect_x1, 0.8, rect_x2, 0.8, lwd = 1.2)
  circos.segments(rect_x1, 0.2, rect_x2, 0.2, lwd = 1.2)
  
  # Left vertical border ONLY on first sector of chromosome
  if (i == fi) {
    circos.segments(CELL_META$xlim[1], 0.2, CELL_META$xlim[1], 0.8, lwd = 1.2)
  }
  
  # Right vertical border ONLY on last sector of chromosome
  if (i == li) {
    circos.segments(CELL_META$xlim[2], 0.2, CELL_META$xlim[2], 0.8, lwd = 1.2)
  }
  
  # Chromosome name only on first sector
  if (i == fi) {
    if (fi != li) {
      center_x <- (CELL_META$xlim[1] + gene_windows$end[li]) / 2
      if (center_x >= CELL_META$xlim[1] && center_x <= CELL_META$xlim[2]) {
        circos.text(center_x, 0.5, cn, cex=0.7, font=2, facing="inside")
      } else {
        circos.text(mean(CELL_META$xlim), 0.5, cn, cex=0.7, font=2, facing="inside")
      }
    } else {
      circos.text(mean(CELL_META$xlim), 0.5, cn, cex=0.7, font=2, facing="inside")
    }
  }
}, track.height = 0.10, bg.border = NA)

circos.clear()
# ---------- 7. LEGENDS ----------
mlv <- intersect(names(mech_colors), unique(gene_meta$Mechanism))
if (length(mlv)==0) mlv <- names(mech_colors)
lg_mech <- Legend(at=mlv, legend_gp=gpar(fill=mech_colors[mlv]),
                  title="Mechanism (gene name)")
lg_pref <- Legend(at=c("Myeloid-biased", "Lymphoid-biased", "Ubiquitous (balanced, no lineage bias)"),
                  legend_gp=gpar(fill=pref_colors),
                  title="Expression preference")
pushViewport(viewport(x=0.90, y=0.5, width=0.18, height=0.85))
lg_list <- list(lg_mech, lg_pref)
for (i in seq_along(lg_list)) {
  pushViewport(viewport(x=0.5, y=1-(i-0.5)/length(lg_list), width=1, height=1/length(lg_list)))
  grid.draw(lg_list[[i]]); popViewport()
}
popViewport()
dev.off()
cat("\n=== DONE ===\nOutput:", output_pdf, "\n")
cat("Per-gene sectors (no blank), within-chr gap=0, chr-block rects overlap 0.5% to kill seams\n")
