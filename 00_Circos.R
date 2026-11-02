# ============================================
# plot_jia_circos(): JIA Figure 1 circos, two layouts, returns stats data.frame
# Input : 0903_Clean case.xlsx (all rows = patients) + HPA scRNA-seq TSVs
# Output: 2 PDFs + gene_meta df (Gene/Mechanism/n_cases/n_mut/novel_class)
# Layout A (outer): chr -> gene -> pref -> case -> mut -> novel -> 8 cell expr
# Layout B (inner): chr -> case -> mut -> novel -> gene -> pref -> 8 cell expr
# ============================================
plot_jia_circos <- function(
    excel_file, hpa_sc_type, hpa_sc_cluster,
    output_pdf_outer = "JIA_Figure1_circos_clinical_outer.pdf",
    output_pdf_inner = "JIA_Figure1_circos_clinical_inner.pdf",
    win_sector = 5e6,
    gene_alias = c("STING1" = "TMEM173"),
    cell_order = c("NK","CD4T","CD8T","Bcell","Monocyte","Macrophage","Neutrophil","DC"),
    mech_colors  = c("deficiency"="#6699CC","gain-of-function"="#FEC260","haploinsufficiency"="#AF478A","na"="grey60"),
    pref_colors  = c("Myeloid-biased"="#82C0AE","Lymphoid-biased"="gray60","Balanced"="#7E75AB"),
    novel_colors = c("Novel"="#BB0021","Reported"="gray80"),
    expr_colors  = c("#007ABA","#FFFFFF","#BB0021"),
    case_colors  = c("#F0F9E8","#7BC47F","#1B7A3D"),
    mut_colors   = c("#FFEFC1","#F49D5C","#F27A49"),
    track_h      = c(chr=0.12, gene=0.12, pref=0.04, case=0.07, mut=0.07, novel=0.05, cell=0.045),
    pdf_w = 16, pdf_h = 12,
    legend_x = 0.88, legend_w = 0.12, legend_h = 0.9) {
  # ---------- packages ----------
  for (pkg in c("openxlsx","circlize","ComplexHeatmap","data.table","dplyr","tidyr","tibble","grid")) {
    if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
      install.packages(pkg, repos="https://cloud.r-project.org"); library(pkg, character.only=TRUE)
    }
  }
  for (pkg in c("org.Hs.eg.db","AnnotationDbi","GenomicFeatures","TxDb.Hsapiens.UCSC.hg38.knownGene")) {
    if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
      if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg); library(pkg, character.only=TRUE)
    }
  }
  select <- dplyr::select
  lymphoid_cells <- cell_order[1:4]; myeloid_cells <- cell_order[5:8]
  # ---------- read excel & per-gene stats ----------
  summary_data <- read.xlsx(excel_file)
  summary_data$gene      <- toupper(trimws(summary_data$gene))
  summary_data$mechanism <- tolower(trimws(summary_data$mechanism))
  summary_data$mechanism[is.na(summary_data$mechanism)] <- "na"
  extract_variants <- function(nucl, hgvs) {
    nucl <- ifelse(is.na(nucl), "", as.character(nucl))
    hgvs <- ifelse(is.na(hgvs), "", as.character(hgvs))
    toks <- trimws(unlist(strsplit(nucl, "/", fixed=TRUE)))
    toks <- toks[toks != "" & toks != "\\"]
    if (length(toks)==0 && hgvs != "") toks <- trimws(unlist(strsplit(hgvs, ":"))[-1])
    toks
  }
  summary_data$mut_set <- lapply(seq_len(nrow(summary_data)), function(i)
    extract_variants(summary_data$Neuclotide[i], summary_data$HGVS[i]))
  summary_data$novel_flag <- tolower(trimws(summary_data$novel))
  gene_meta <- summary_data %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(Mechanism = paste(unique(mechanism[mechanism != ""]), collapse="/"),
                     n_cases = n(),
                     n_mut = length(unique(unlist(mut_set))),
                     novel_class = ifelse(any(grepl("novel", novel_flag)), "Novel", "Reported"),
                     .groups="drop") %>%
    dplyr::mutate(Gene=as.character(gene),
                  Mechanism=ifelse(Mechanism==""|is.na(Mechanism),"na",Mechanism)) %>%
    dplyr::select(Gene, Mechanism, n_cases, n_mut, novel_class)
  gene_meta$line_color  <- unname(mech_colors[gene_meta$Mechanism])
  gene_meta$line_color[is.na(gene_meta$line_color)] <- "grey60"
  gene_meta$novel_color <- novel_colors[gene_meta$novel_class]
  cat("\n=== Per-gene stats (returned df) ===\n")
  print(data.frame(Gene=gene_meta$Gene, Mechanism=gene_meta$Mechanism,
                   n_cases=gene_meta$n_cases, n_mut=gene_meta$n_mut,
                   novel=gene_meta$novel_class))
  # ---------- hg38 coordinates ----------
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
  # ---------- HPA single-cell expression ----------
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
  # ---------- expression preference ----------
  lymphoid_mean <- rowMeans(mat_sc_log[, lymphoid_cells])
  myeloid_mean  <- rowMeans(mat_sc_log[, myeloid_cells])
  pref_score <- (myeloid_mean - lymphoid_mean) / (myeloid_mean + lymphoid_mean)
  pref_score[is.nan(pref_score)] <- 0
  pref_class <- ifelse(pref_score > 0.15, "Myeloid-biased",
                       ifelse(pref_score < -0.15, "Lymphoid-biased", "Balanced"))
  gene_meta$pref_score <- pref_score
  gene_meta$pref_class <- pref_class
  gene_meta$pref_color <- pref_colors[pref_class]
  # ---------- sectors & color functions ----------
  gene_windows <- data.frame(chr=gene_meta$chr,
                             start=pmax(gene_meta$start - win_sector, 0),
                             end=gene_meta$end + win_sector,
                             stringsAsFactors=FALSE)
  gaps <- c()
  for (i in 1:(nrow(gene_meta)-1)) {
    if (gene_meta$chr[i] == gene_meta$chr[i+1]) gaps <- c(gaps, 0) else gaps <- c(gaps, 6)
  }
  gaps <- c(gaps, 8)
  chr_first <- c(); chr_last <- c()
  for (cn in unique(gene_meta$chr)) {
    idx <- which(gene_meta$chr == cn)
    chr_first[cn] <- min(idx); chr_last[cn] <- max(idx)
  }
  max_sc <- max(mat_sc_log, na.rm=TRUE); if (max_sc <= 0) max_sc <- 1
  expr_col_fun <- colorRamp2(breaks=c(0, max_sc*0.5, max_sc), colors=expr_colors)
  max_cases <- max(gene_meta$n_cases); if (max_cases < 2) max_cases <- 2
  max_mut   <- max(gene_meta$n_mut);   if (max_mut   < 2) max_mut   <- 2
  case_col_fun <- colorRamp2(breaks=c(1, (1+max_cases)/2, max_cases), colors=case_colors)
  mut_col_fun  <- colorRamp2(breaks=c(1, (1+max_mut)/2, max_mut),   colors=mut_colors)
  th <- as.list(track_h)
  # ---------- track drawing helpers ----------
  tr_cell <- function(ct) {
    local({
      ct_local <- ct
      circos.track(ylim=c(0,1), panel.fun=function(x,y) {
        i <- CELL_META$sector.numeric.index
        v <- as.numeric(mat_sc_log[i, ct_local])
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                    col=expr_col_fun(v), border="black", lwd=0.5)
        if (i == 1) circos.text(CELL_META$xlim[1], 0.5, ct_local, cex=0.55, font=1,
                                facing="inside", adj=c(0,0.5), col="black")
      }, track.height=th$cell, bg.border=NA)
    })
  }
  tr_novel <- function() {
    circos.track(ylim=c(0,1), panel.fun=function(x,y) {
      i <- CELL_META$sector.numeric.index
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col=gene_meta$novel_color[i], border="black", lwd=0.5)
    }, track.height=th$novel, bg.border=NA)
  }
  tr_mut <- function() {
    circos.track(ylim=c(0,1), panel.fun=function(x,y) {
      i <- CELL_META$sector.numeric.index
      v <- gene_meta$n_mut[i]
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col=mut_col_fun(v), border="black", lwd=0.5)
      circos.text(mean(CELL_META$xlim), 0.5, as.character(v), cex=0.6, font=2,
                  facing="inside", col=ifelse(v > (1+max_mut)/2, "white", "black"))
    }, track.height=th$mut, bg.border=NA)
  }
  tr_case <- function() {
    circos.track(ylim=c(0,1), panel.fun=function(x,y) {
      i <- CELL_META$sector.numeric.index
      v <- gene_meta$n_cases[i]
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col=case_col_fun(v), border="black", lwd=0.5)
      circos.text(mean(CELL_META$xlim), 0.5, as.character(v), cex=0.6, font=2,
                  facing="inside", col=ifelse(v > (1+max_cases)/2, "white", "black"))
    }, track.height=th$case, bg.border=NA)
  }
  tr_pref <- function() {
    circos.track(ylim=c(0,1), panel.fun=function(x,y) {
      i <- CELL_META$sector.numeric.index
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col=gene_meta$pref_color[i], border="black", lwd=0.5)
      short <- ifelse(gene_meta$pref_class[i]=="Myeloid-biased", "M",
                      ifelse(gene_meta$pref_class[i]=="Lymphoid-biased", "L", "B"))
      circos.text(mean(CELL_META$xlim), 0.5, short, cex=0.55, font=2,
                  facing="inside", col="white")
    }, track.height=th$pref, bg.border=NA)
  }
  tr_gene <- function() {
    circos.track(ylim=c(0,1), panel.fun=function(x,y) {
      i <- CELL_META$sector.numeric.index
      circos.text(mean(CELL_META$xlim), 0.5, gene_meta$Gene[i],
                  cex=0.9, font=1, facing="inside", col=gene_meta$line_color[i])
    }, track.height=th$gene, bg.border=NA)
  }
  tr_chr <- function() {
    circos.track(ylim=c(0,1), panel.fun=function(x,y) {
      i <- CELL_META$sector.numeric.index
      cn <- gene_meta$chr[i]; fi <- chr_first[cn]; li <- chr_last[cn]
      w <- CELL_META$xlim[2] - CELL_META$xlim[1]
      ov <- w * 0.005
      rect_x1 <- CELL_META$xlim[1] - ov
      rect_x2 <- CELL_META$xlim[2] + ov
      circos.rect(rect_x1, 0.2, rect_x2, 0.8, col="grey85", border=NA)
      circos.segments(rect_x1, 0.8, rect_x2, 0.8, lwd=1.2)
      circos.segments(rect_x1, 0.2, rect_x2, 0.2, lwd=1.2)
      if (i == fi) circos.segments(CELL_META$xlim[1], 0.2, CELL_META$xlim[1], 0.8, lwd=1.2)
      if (i == li) circos.segments(CELL_META$xlim[2], 0.2, CELL_META$xlim[2], 0.8, lwd=1.2)
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
    }, track.height=th$chr, bg.border=NA)
  }
  # ---------- single-plot driver ----------
  draw_one <- function(output_pdf, clinical_pos) {
    circos.clear()
    pdf(output_pdf, width=pdf_w, height=pdf_h)
    circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0), gap.after=gaps)
    fa <- paste0(gene_windows$chr, ":", gene_windows$start, "-", gene_windows$end)
    circos.initialize(factors=factor(fa, levels=fa),
                      xlim=cbind(gene_windows$start, gene_windows$end))
    if (clinical_pos == "outer") {
      for (ct in rev(cell_order)) tr_cell(ct)
      tr_novel(); tr_mut(); tr_case(); tr_pref(); tr_gene(); tr_chr()
    } else {
      for (ct in rev(cell_order)) tr_cell(ct)
      tr_pref(); tr_gene(); tr_novel(); tr_mut(); tr_case(); tr_chr()
    }
    mlv <- intersect(names(mech_colors), unique(gene_meta$Mechanism))
    if (length(mlv)==0) mlv <- names(mech_colors)
    lg_mech  <- Legend(at=mlv, legend_gp=gpar(fill=mech_colors[mlv]), title="Mechanism (gene name)",
                       title_gp=gpar(fontsize=9, fontface="bold"), labels_gp=gpar(fontsize=8))
    lg_pref  <- Legend(at=names(pref_colors), legend_gp=gpar(fill=pref_colors), title="Expression preference",
                       title_gp=gpar(fontsize=9, fontface="bold"), labels_gp=gpar(fontsize=8))
    lg_sc    <- Legend(at=c(0, round(max_sc/2,1), round(max_sc,1)), col_fun=expr_col_fun,
                       title="scExpression (log2 nTPM+1)",
                       title_gp=gpar(fontsize=9, fontface="bold"), labels_gp=gpar(fontsize=8))
    lg_case  <- Legend(at=c(1, round((1+max_cases)/2,1), max_cases), col_fun=case_col_fun,
                       title="Cases per gene",
                       title_gp=gpar(fontsize=9, fontface="bold"), labels_gp=gpar(fontsize=8))
    lg_mut   <- Legend(at=c(1, round((1+max_mut)/2,1), max_mut), col_fun=mut_col_fun,
                       title="Mutations per gene",
                       title_gp=gpar(fontsize=9, fontface="bold"), labels_gp=gpar(fontsize=8))
    lg_novel <- Legend(at=names(novel_colors), legend_gp=gpar(fill=novel_colors), title="Novel variant",
                       title_gp=gpar(fontsize=9, fontface="bold"), labels_gp=gpar(fontsize=8))
    pushViewport(viewport(x=legend_x, y=0.5, width=legend_w, height=legend_h))
    lg_list <- list(lg_mech, lg_pref, lg_sc, lg_case, lg_mut, lg_novel)
    for (i in seq_along(lg_list)) {
      pushViewport(viewport(x=0.5, y=1-(i-0.5)/length(lg_list), width=1, height=1/length(lg_list)))
      grid.draw(lg_list[[i]]); popViewport()
    }
    popViewport()
    dev.off()
    cat("PDF:", output_pdf, "\n")
  }
  # ---------- render both layouts ----------
  draw_one(output_pdf_outer, "outer")
  draw_one(output_pdf_inner, "inner")
  cat("\n=== DONE ===\n")
  gene_meta
}


res <- plot_jia_circos(
  excel_file     = "E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx",
  hpa_sc_type    = "E:\\Cohort PPT\\JIA\\rna_single_cell_type.tsv",
  hpa_sc_cluster = "E:\\Cohort PPT\\JIA\\rna_single_cell_cluster.tsv"
)
print(res)


###check again
library(openxlsx)
metainfo<-read.xlsx("E:\\Cohort PPT\\JIA\\0903_Clean case.xlsx")
table(metainfo$gene)

genes<-unique(metainfo$gene)
check_mutclini<-data.frame()
for(i in 1:length(genes)){
  data1<-metainfo[metainfo$gene==genes[i],] 
  mut<-unique(as.character(unlist(sapply(data1$Protein,function(x){unlist(strsplit(x,split = "/"))}))))
  check_mutclini<-rbind(
    check_mutclini,
    data.frame(
      gene=genes[i],
      n_case=nrow(data1),
      n_mut=length(mut))
  )
}


