library(Seurat)
library(patchwork)
library(ggplot2)
library(AUCell)
library(ggsignif)
library(ggpubr)
library(ggplotify)
library(msigdbr)
library(dplyr)

setwd("E:/Cohort PPT/JIA/code")
pbmc1<-readRDS("pbmc1.RDS")

#####1 Differentially expressed genes across cell types, grouped by different gene mutations
Mut_DEG <- function(seu_obj, mut_col = "mutation",
                    refs = c("Healthy controls","Polygenic JIA"),
                    celltypes = NULL, exclude_celltypes = c("Erythrocyte","Platelet"),
                    min_cells = 10, min_cells_grp = 3, min_pct = 0.1, logfc = 0, verbose = TRUE){
  library(Seurat); library(dplyr); library(tibble)
  md <- seu_obj@meta.data
  all_ct <- unique(as.character(md$celltype))
  if(is.null(celltypes)) celltypes <- setdiff(all_ct, exclude_celltypes)
  celltypes <- intersect(celltypes, all_ct)
  ref_tag <- c("Healthy controls" = "HC", "Polygenic JIA" = "Poly")
  all_mut <- setdiff(unique(as.character(md[[mut_col]])), refs)
  total <- length(celltypes) * length(all_mut) * length(refs)
  out <- list(); skip <- list(); k <- 0; nrow_acc <- 0
  t0 <- Sys.time()
  if(verbose) cat(sprintf("PLAN: %d celltypes x %d mutations x %d refs = %d tests\n",
                          length(celltypes), length(all_mut), length(refs), total), file = stderr())
  for(ct in celltypes){
    idx <- which(as.character(md$celltype) == ct)
    muts <- setdiff(unique(as.character(md[[mut_col]][idx])), refs)
    if(length(idx) < min_cells){
      k <- k + length(muts) * length(refs)
      skip[[length(skip)+1]] <- data.frame(celltype = ct, mutation = NA, ref = NA, n = length(idx), reason = "celltype cells < min_cells")
      if(verbose) cat(sprintf("SKIP  [%3d/%3d] %-18s | only %d cells (< %d)\n", k, total, ct, length(idx), min_cells), file = stderr())
      next
    }
    sub <- seu_obj[, rownames(md)[idx]]
    Idents(sub) <- mut_col
    for(m in muts){
      for(r in refs){
        k <- k + 1
        n1 <- sum(as.character(sub@meta.data[[mut_col]]) == m)
        n2 <- sum(as.character(sub@meta.data[[mut_col]]) == r)
        if(n1 < min_cells_grp || n2 < min_cells_grp){
          skip[[length(skip)+1]] <- data.frame(celltype = ct, mutation = m, ref = r, n = paste0(n1, " vs ", n2), reason = "group cells < min_cells_grp")
          if(verbose) cat(sprintf("SKIP  [%3d/%3d] %-18s | %-22s vs %-4s | n = %d vs %d\n", k, total, ct, m, ref_tag[r], n1, n2), file = stderr())
          next
        }
        if(verbose) cat(sprintf("RUN   [%3d/%3d] %-18s | %-22s vs %-4s | n = %d vs %d\n", k, total, ct, m, ref_tag[r], n1, n2), file = stderr())
        t1 <- Sys.time()
        mk <- FindMarkers(sub, ident.1 = m, ident.2 = r, only.pos = FALSE,
                          min.pct = min_pct, logfc.threshold = logfc, verbose = verbose)
        dt <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
        if(nrow(mk) == 0){
          if(verbose) cat(sprintf("EMPTY [%3d/%3d] %-18s | %-22s vs %-4s\n", k, total, ct, m, ref_tag[r]), file = stderr())
          next
        }
        mk <- mk %>% rownames_to_column("gene") %>%
          mutate(celltype = ct, mutation = m, ref = r,
                 comparison = paste0(m, "_vs_", ref_tag[r]),
                 n_mut = n1, n_ref = n2)
        out[[length(out)+1]] <- mk
        nrow_acc <- nrow_acc + nrow(mk)
        el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        eta <- el / k * (total - k)
        if(verbose) cat(sprintf("DONE  [%3d/%3d] %-18s | %-22s vs %-4s | %5d genes | %5.1fs | rows %6d | ETA %4.1f min\n",
                                k, total, ct, m, ref_tag[r], nrow(mk), dt, nrow_acc, eta/60), file = stderr())
      }
    }
  }
  el_all <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if(verbose) cat(sprintf("\nEND: %d combos | %d rows | %d skipped | %.1f min\n", k, nrow_acc, length(skip), el_all/60), file = stderr())
  list(deg = if(length(out) > 0) bind_rows(out) else NULL,
       skipped = if(length(skip) > 0) bind_rows(skip) else NULL)
}
cells<-cells[-match(c("Platelet","Erythrocyte"),cells)]
res_mut <- Mut_DEG(pbmc1, celltypes =cells)
saveRDS(res_mut,file="1109_DEG.RDS")

DEGs<-readRDS("E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\1109_DEG.RDS")

####(1) DEG count for each celltype, VS healthy control, VS polygenic JIA patients
library(dplyr); library(tidyr); library(ggplot2)
deg_sig <- DEGs$deg %>%
  filter(p_val_adj < 0.01, abs(avg_log2FC) > 1.5) %>%
  mutate(dir = ifelse(avg_log2FC > 0, "up", "down"))

cnt <- deg_sig %>% count(celltype, mutation, ref, dir) %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0)
if(!"up"   %in% names(cnt)) cnt$up   <- 0
if(!"down" %in% names(cnt)) cnt$down <- 0

cnt <- expand_grid(celltype = unique(DEGs$deg$celltype),
                   mutation = unique(DEGs$deg$mutation),
                   ref      = unique(DEGs$deg$ref)) %>%
  left_join(cnt, by = c("celltype", "mutation", "ref")) %>%
  mutate(up = replace_na(up, 0), down = replace_na(down, 0), total = up + down) %>%
  arrange(ref, mutation, celltype)

CT_ALL <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","Cytotoxic CD8 T","NK","Treg",
            "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC",
            "Platelet","Erythrocyte")
CT_COL_ALL <- c("#847AB3","#96AF95","#6699CC","#9CD2ED","#86C7B4","#007ABA",
                "#F1BAAB","#BB0021","#D9C2D9","#FEC260","#AF478A",
                "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")
names(CT_COL_ALL) <- CT_ALL
CT_USE <- intersect(CT_ALL, unique(DEGs$deg$celltype))
CT_COL <- CT_COL_ALL[CT_USE]

MUTS <- unique(DEGs$deg$mutation)
MUT_LAB <- setNames(vapply(MUTS, function(m){
  if(grepl("366", m)) "LACC1"
  else if(grepl("781", m)) "LPIN2"
  else if(grepl("513", m)) "NOD2 (M513K)"
  else if(grepl("495", m)) "NOD2 (C495Y)"
  else if(grepl("250", m)) "PSTPIP1"
  else m
}, character(1)), MUTS)
MUT_LEVELS <- c("LACC1","LPIN2","NOD2 (M513K)","NOD2 (C495Y)","PSTPIP1")

plot_ref <- function(df, ref_name, mode = c("dir","celltype"),
                     col_up = "#EC706E", col_down = "#234091",
                     ct_colors = NULL, light_alpha = 0.35, bar_w = 0.7,
                     facet_nrow = 1){
  mode <- match.arg(mode)
  if(is.null(ct_colors)) ct_colors <- CT_COL
  d <- df %>% filter(ref == ref_name) %>%
    pivot_longer(c(up, down), names_to = "dir", values_to = "n") %>%
    mutate(val = ifelse(dir == "up", n, -n),
           celltype = factor(celltype, levels = CT_USE),
           dir = factor(dir, levels = c("up","down")))
  d$mut_show <- factor(unname(MUT_LAB[d$mutation]), levels = MUT_LEVELS)
  if(mode == "dir"){
    p <- ggplot(d, aes(x = celltype, y = val, fill = dir)) +
      geom_col(width = bar_w) +
      scale_fill_manual(values = c(up = col_up, down = col_down), labels = c("Up","Down"))
  } else {
    p <- ggplot(d, aes(x = celltype, y = val)) +
      geom_col(data = d %>% filter(dir == "up"), aes(fill = celltype), width = bar_w) +
      geom_col(data = d %>% filter(dir == "down"), aes(fill = celltype), width = bar_w, alpha = light_alpha) +
      scale_fill_manual(values = ct_colors)
  }
  p + geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
    facet_wrap(~ mut_show, nrow = facet_nrow) +
    scale_y_continuous(labels = abs) +
    labs(x = NULL, y = "DEG count", fill = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, colour = "black"),
          axis.text.y = element_text(size = 7, colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 7.5, colour = "black"),
          panel.border = element_rect(colour = "gray90", fill = NA, linetype = "dashed", linewidth = 0.3),
          panel.spacing.x = unit(0.35, "lines"),
          legend.position = "bottom")
}

p_hc_a   <- plot_ref(cnt, "Healthy controls",  mode = "dir")
p_poly_a <- plot_ref(cnt, "Polygenic JIA",      mode = "dir")
p_hc_b   <- plot_ref(cnt, "Healthy controls",  mode = "celltype")
p_poly_b <- plot_ref(cnt, "Polygenic JIA",      mode = "celltype")
p_hc_a; p_poly_a; p_hc_b; p_poly_b

deg_sig$gene_lab <- unname(MUT_LAB[deg_sig$mutation])

chk_gene <- deg_sig %>%
  group_by(gene_lab, ref, dir) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0)
if(!"up"   %in% names(chk_gene)) chk_gene$up   <- 0
if(!"down" %in% names(chk_gene)) chk_gene$down <- 0
chk_gene <- chk_gene %>% mutate(total = up + down) %>% arrange(ref, gene_lab)
chk_from_cnt <- cnt %>%
  mutate(gene_lab = unname(MUT_LAB[mutation])) %>%
  group_by(gene_lab, ref) %>%
  summarise(up = sum(up), down = sum(down), total = sum(total), .groups = "drop") %>%
  arrange(ref, gene_lab)

chk_gene
chk_from_cnt

library(patchwork)
p_hc_b/ p_poly_b


####(2) Check the intersection between DEGS and signature genes
library(dplyr); library(tidyr); library(tibble); library(rlang); library(msigdbr)
load("E:\\通路基因集合\\gene4pathway.Rdata")
IFN28 <- c("DDX60","EPSTI1","HERC5","HERC6","IFI27","IFI44","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFIT5","ISG15","LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SIGLEC1","SPATS2L","USP18","CXCL10","GBP1","SOCS1")
pyrop <- c("AIM2","APIP","CASP1","CASP4","CASP8","DHX9","ELANE","GSDMA","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB","NAIP","NLRC4","NLRP1","NLRP6","NLRP9","TREM2","ZBP1","NLRP3","PYCARD","NEK7","IL1B")
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
auto_genes <- reactome %>% filter(gs_name == "REACTOME_AUTOPHAGY") %>% pull(gene_symbol) %>% unique()
hm <- msigdbr(species = "Homo sapiens", category = "H")
gly_genes <- hm %>% filter(gs_name == "HALLMARK_GLYCOLYSIS") %>% pull(gene_symbol) %>% unique()
oxphos_genes <- hm %>% filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% pull(gene_symbol) %>% unique()
arachidonic_genes <- c("PLA2G4A","PLA2G6","PLA2G2A","PTGS1","PTGS2","ALOX5","ALOX5AP","ALOX12","ALOX15","ALOX15B","LTC4S","PTGES","PTGES2","TBXAS1","GPX4","CYP2J2","CYP4F2","EPHX2","HPGDS")
inflam_genes <- hm %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(gene_symbol) %>% unique()

gs <- list(IFN28 = IFN28, pyrop = pyrop, autophagy = auto_genes, glycolysis = gly_genes,
           OXPHOS = oxphos_genes, arachidonic = arachidonic_genes, inflam = inflam_genes,
           NFKB = NFKB$gene)

hit_ct_dir <- deg_sig %>%
  group_by(ref, celltype, dir) %>%
  group_modify(~ tibble(n_deg = nrow(.x), !!!lapply(gs, function(s) sum(.x$gene %in% s)))) %>%
  ungroup()

hit_hc   <- hit_ct_dir %>% filter(ref == "Healthy controls") %>% arrange(celltype, dir)
hit_poly <- hit_ct_dir %>% filter(ref == "Polygenic JIA")     %>% arrange(celltype, dir)
hit_hc
hit_poly


#####2 GSEA analysis based on DEGs derived from above results, VS HC; VS polygenic JIA patients for each cell type
####(1) GSEA running
library(dplyr); library(tidyr); library(tibble); library(clusterProfiler); library(msigdbr)
fix_col <- function(x){ if(!"gene_symbol" %in% names(x)) names(x)[names(x) == "gene"] <- "gene_symbol"; x }
hall <- fix_col(msigdbr(species = "Homo sapiens", category = "H"))
kegg <- fix_col(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"))
if(nrow(kegg) == 0) kegg <- fix_col(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_MEDICUS"))

t2g_list <- list(HALLMARK = hall[, c("gs_name", "gene_symbol")],
                 KEGG     = kegg[, c("gs_name", "gene_symbol")])

run_gsea <- function(deg_all, t2g_list, pvalueCutoff = 1, minGSSize = 10,
                     maxGSSize = 1000, nPermSimple = 1000){
  comps <- deg_all %>% distinct(celltype, mutation, ref, comparison) %>% arrange(comparison, celltype)
  objs <- list(); dfs <- list(); t0 <- Sys.time(); n <- nrow(comps)
  for(i in seq_len(n)){
    ct <- comps$celltype[i]; mu <- comps$mutation[i]
    rf <- comps$ref[i];      cp <- comps$comparison[i]
    sub <- deg_all %>% filter(celltype == ct, mutation == mu, ref == rf) %>%
      filter(!is.na(avg_log2FC)) %>% arrange(desc(avg_log2FC))
    if(nrow(sub) < 10) next
    rk <- sub$avg_log2FC
    names(rk) <- sub$gene
    rk <- rk[!duplicated(names(rk))]
    rk <- sort(rk, decreasing = TRUE)
    for(db in names(t2g_list)){
      key <- paste(ct, cp, db, sep = "|")
      g <- GSEA(rk, TERM2GENE = t2g_list[[db]], pvalueCutoff = pvalueCutoff,
                minGSSize = minGSSize, maxGSSize = maxGSSize,
                nPermSimple = nPermSimple, verbose = FALSE)
      if(is.null(g) || nrow(g@result) == 0) next
      objs[[key]] <- g
      df <- as.data.frame(g@result)
      df$celltype <- ct; df$mutation <- mu; df$ref <- rf
      df$comparison <- cp; df$database <- db
      dfs[[key]] <- df
      cat(sprintf("DONE [%3d/%3d] %-16s %-22s %-8s | %4d terms\n",
                  i, n, ct, cp, db, nrow(df)), file = stderr())
    }
  }
  gsea_res <- bind_rows(dfs)
  list(gsea_obj = objs, gsea_res = gsea_res)
}

gsea_all <- run_gsea(DEGs$deg, t2g_list)
gsea_res <- gsea_all$gsea_res
gsea_obj <- gsea_all$gsea_obj
save(gsea_all,gsea_res,gsea_obj,file="E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\1111_GSEA.Rdata")
openxlsx::write.xlsx(gsea_res,file="1111_MutgroupGSEA.xlsx")

####(2) bubble plot for GSEA results
MUTS <- unique(gsea_res$mutation)
MUT_LAB <- setNames(vapply(MUTS, function(m){
  if(grepl("366", m)) "LACC1" else if(grepl("781", m)) "LPIN2"
  else if(grepl("513", m)) "NOD2 (M513K)" else if(grepl("495", m)) "NOD2 (C495Y)"
  else if(grepl("250", m)) "PSTPIP1" else m
}, character(1)), MUTS)
MUT_LEVELS <- c("LACC1", "LPIN2", "NOD2 (M513K)", "NOD2 (C495Y)", "PSTPIP1")

plot_gsea_bubble <- function(data, pathway_order, pathway_labels, ct_use,
                             mode = c("nes_size", "padj_size"),
                             fill_palette = col2,
                             shape_type = 16, size_range = c(2, 6),
                             size_limits = NULL, angle = 45,
                             padj_cut = 0.01, nes_cut = 1.5,
                             ref_levels = c("Healthy controls", "Polygenic JIA"),
                             ref_labels = c("vs HC", "vs Poly"),
                             ref_split = TRUE, facet_drop = FALSE){
  library(ggplot2); library(dplyr)
  mode <- match.arg(mode)
  d <- data %>%
    filter(ID %in% pathway_order, celltype %in% ct_use) %>%
    mutate(NES = as.numeric(NES), padj = as.numeric(p.adjust),
           neglog10padj = -log10(padj + 1e-300)) %>%
    filter(padj < padj_cut, NES > nes_cut) %>%
    mutate(ID = factor(ID, levels = pathway_order, labels = pathway_labels),
           celltype = factor(celltype, levels = ct_use),
           ref = factor(ref, levels = ref_levels, labels = ref_labels),
           mut = factor(unname(MUT_LAB[mutation]), levels = MUT_LEVELS))
  size_var <- if(mode == "nes_size") "NES" else "neglog10padj"
  d$sz <- d[[size_var]]
  if(!is.null(size_limits)) d$sz <- pmin(pmax(d$sz, size_limits[1]), size_limits[2])
  base_theme <- theme_bw() +
    theme(axis.text = element_text(color = "black", size = 8),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(color = "black"),
          strip.text = element_text(color = "black", size = 9),
          legend.title = element_text(color = "black", size = 8),
          legend.text = element_text(color = "black", size = 7),
          panel.grid.major = element_line(color = "grey85", linewidth = 0.3, linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.5))
  mk <- function(dd, xvar, yvar, facet){
    dd <- droplevels(dd)
    if(mode == "nes_size"){
      p <- ggplot(dd, aes(x = .data[[xvar]], y = .data[[yvar]],
                          size = sz, colour = neglog10padj)) +
        scale_size_continuous(range = size_range) +
        scale_colour_gradientn(colors = fill_palette) +
        labs(x = "", y = "", size = "NES", colour = "-log10(p.adjust)")
    } else {
      p <- ggplot(dd, aes(x = .data[[xvar]], y = .data[[yvar]],
                          size = sz, colour = NES)) +
        scale_size_continuous(range = size_range) +
        scale_colour_gradientn(colors = fill_palette) +
        labs(x = "", y = "", size = "-log10(p.adjust)", colour = "NES")
    }
    p + geom_point(shape = shape_type, stroke = 0.5) + facet + base_theme +
      theme(axis.text.x = element_text(angle = angle, hjust = 1))
  }
  pg <- list(mut_row = mk(d, "ID", "celltype", facet_grid(mut ~ ref, drop = facet_drop)),
             mut_col = mk(d, "celltype", "ID", facet_grid(ref ~ mut, drop = facet_drop)))
  sp <- NULL
  if(ref_split){
    sp <- lapply(ref_labels, function(r){
      dd <- d[d$ref == r, , drop = FALSE]
      list(mut_row = mk(dd, "ID", "celltype", facet_grid(mut ~ ., drop = TRUE)),
           mut_col = mk(dd, "celltype", "ID", facet_grid(. ~ mut, drop = TRUE)))
    })
    names(sp) <- ref_labels
  }
  list(split = sp, combined = pg, data = d)
}

col2 <- colorRampPalette(c("#007ABA", "#F7F7F7", "#BB0021"))(101)

pathway_order <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING")
pathway_labels <- c("TNF via NF-κB", "Inflammatory response",
                    "Interferon-γ response", "Interferon-α response",
                    "MAPK signaling", "T cell receptor signaling", "IL6-JAK-STAT3")
ct_use <- c("CD14 Monocyte", "CD16 Monocyte", "pDC", "Naive CD4 T", "Cytotoxic CD8 T")

p1 <- plot_gsea_bubble(
  data = gsea_res, pathway_order = pathway_order, pathway_labels = pathway_labels,
  ct_use = ct_use, mode = "nes_size",
  fill_palette = col2,
  shape_type = 15, size_range = c(2.5, 5.5),
  size_limits = c(1.5, 2.8), angle = 45,
  padj_cut = 0.01, nes_cut = 1.5,
  ref_levels = c("Healthy controls", "Polygenic JIA"), ref_labels = c("vs HC", "vs Poly"),
  ref_split = TRUE, facet_drop = FALSE)
p1$split[["vs HC"]]$mut_row
p1$split[["vs Poly"]]$mut_row
p1$split[["vs HC"]]$mut_col
p1$split[["vs Poly"]]$mut_col

p2 <- plot_gsea_bubble(
  data = gsea_res, pathway_order = pathway_order, pathway_labels = pathway_labels,
  ct_use = ct_use, mode = "padj_size",
  fill_palette = col2,
  shape_type = 16, size_range = c(2.5, 5.5),
  size_limits = c(1.3, 4), angle = 45,
  padj_cut = 0.01, nes_cut = 1.5,
  ref_levels = c("Healthy controls", "Polygenic JIA"), ref_labels = c("vs HC", "vs Poly"),
  ref_split = TRUE, facet_drop = FALSE)
p2$split[["vs HC"]]$mut_row
p2$split[["vs Poly"]]$mut_row
