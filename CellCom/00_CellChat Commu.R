library(Seurat); library(CellChat); library(dplyr); library(patchwork); library(circlize); library(future)
setwd("E:\\Cohort PPT\\JIA\\code\\CellCom")


pbmc1<-readRDS("E:\\Cohort PPT\\JIA\\code\\pbmc1.RDS")
NBoot <- 30; MinCells <- 10; n_workers <- 8

ct_keep <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
             "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC")
color_ct <- c("Naive CD4 T"="#847AB3","Memory CD4 T"="#96AF95","Naive CD8 T"="#6699CC",
              "NK"="#86C7B4","Treg"="#007ABA","Cytotoxic CD8 T"="#9CD2ED","γδT"="#F1BAAB",
              "Naive B"="#BB0021","Memory B"="#D9C2D9","Plasma"="#FEC260",
              "CD14 Monocyte"="#AF478A","CD16 Monocyte"="#EC706E","pDC"="#B5BBE3")
mut_levels <- c("Healthy controls","Polygenic JIA","S366Ter/N19KfsTer31",
                "I67R/K781TfsTer7","M513K","C495Y","E250K")
short_lab <- c("HC","Poly","S366Ter","I67R","M513K","C495Y","E250K"); names(short_lab) <- mut_levels
color_grp <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#B07AA1","#EDC948"); names(color_grp) <- short_lab

######1 Proprecessing and group
pbmc1 <- subset(pbmc1, celltype %in% ct_keep)
gc();unique(pbmc1$celltype)
pbmc1$celltype <- factor(pbmc1$celltype, levels = ct_keep)
pbmc1$group_chat <- factor(pbmc1$mutation, levels = mut_levels)
pbmc1 <- pbmc1[, which(!is.na(pbmc1$group_chat))]

print(table(pbmc1$datasets, pbmc1$mutation))
n_tab <- table(pbmc1$celltype, pbmc1$group_chat)
print(addmargins(n_tab))
print(apply(n_tab, 2, function(x) sum(x >= MinCells)))

######3 Cellchat subject construction
options(future.globals.maxSize = 4 * 1024^3)
setwd("E:\\Cohort PPT\\JIA\\code\\CellCom")
chat_list <- list()
for (g in mut_levels) {
  message("[", format(Sys.time(), "%H:%M:%S"), "] START ", short_lab[g])
  cells_g <- colnames(pbmc1)[which(pbmc1$group_chat == g)]
  if (!length(cells_g)) next
  sub <- pbmc1[, cells_g]
  if (length(unique(sub$celltype)) < 2) next
  print(table(sub$celltype))
  cc <- createCellChat(GetAssayData(sub, assay = "RNA", slot = "data"), meta = sub@meta.data, group.by = "celltype")
  cc@idents <- droplevels(cc@idents)
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- projectData(cc, PPI.human)
  cc <- computeCommunProb(cc, raw.use = FALSE, population.size = TRUE, nboot = NBoot)
  cc <- filterCommunication(cc, min.cells = MinCells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  chat_list[[short_lab[g]]] <- cc
  saveRDS(cc, paste0("chat_", short_lab[g], ".RDS"))
  message("[", format(Sys.time(), "%H:%M:%S"), "] DONE  ", short_lab[g], " | cells: ", ncol(sub))
  rm(sub); gc()
}
print(names(chat_list))

######4 merge
chat_list <- lapply(chat_list, liftCellChat, group.new = ct_keep)
chat_merge <- mergeCellChat(chat_list, add.names = names(chat_list), cell.prefix = TRUE)
class(chat_merge)
saveRDS(chat_merge, "JIA_cellchat_mutation7_merge.RDS")
names(chat_merge@net)
str(chat_merge@net, max.level = 1)
slotNames(chat_merge)
names(chat_list$HC@net)
dim(chat_merge@net$HC$count)
packageVersion("CellChat")

######5 Overall interaction number and strength
col_one <- "#6699CC"
grp <- names(chat_merge@net)
df <- data.frame(
  group  = factor(grp, levels = grp),
  count  = sapply(grp, function(g) sum(chat_merge@net[[g]]$count)),
  weight = sapply(grp, function(g) sum(chat_merge@net[[g]]$weight))
)

hc_i   <- which(grepl("^HC", grp, ignore.case = TRUE))
poly_i <- which(grepl("poly", grp, ignore.case = TRUE))

make_hline <- function(y_hc, y_poly){
  c(lapply(y_hc,   function(y) geom_hline(yintercept = y, linetype = "dashed", colour = "grey30",  linewidth = 0.3)),
    lapply(y_poly, function(y) geom_hline(yintercept = y, linetype = "dashed", colour = "#EC706E", linewidth = 0.3)))
}

p1 <- ggplot(df, aes(x = group, y = count)) +
  geom_col(width = 0.65, fill = col_one, colour = "black", linewidth = 0.3) +
  make_hline(df$count[hc_i], df$count[poly_i]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  labs(x = NULL, y = "Number of interactions")

p2 <- ggplot(df, aes(x = group, y = weight)) +
  geom_col(width = 0.65, fill = col_one, colour = "black", linewidth = 0.3) +
  make_hline(df$weight[hc_i], df$weight[poly_i]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  labs(x = NULL, y = "Interaction strength")

p1 | p2

######6 rankNet, fingding different pathway across groups
length(chat_merge@netP$HC$pathways)
length(chat_merge@netP$Poly$pathways) #different pathway numbers among different samples

####(1) primary version
rankNet(chat_merge, mode = "comparison", measure = "weight", comparison = seq_along(chat_list), stacked = TRUE, do.stat = FALSE)
rankNet(chat_merge, mode = "comparison", measure = "count", comparison = seq_along(chat_list), stacked = TRUE, do.stat = FALSE)
rankNet(chat_merge, mode = "comparison", measure = "weight", comparison = seq_along(chat_list), stacked = FALSE, do.stat = TRUE)

####(2) extract intersected pathways
pw <- Reduce(intersect, lapply(names(chat_merge@netP), 
function(g) chat_merge@netP[[g]]$pathways))
length(pw)
flow <- sapply(names(chat_merge@netP), function(g) {
  x <- chat_merge@netP[[g]]
  apply(x$prob, 3, sum)[match(pw, x$pathways)]
})
rownames(flow) <- pw
flow_rel <- sweep(flow, 1, rowSums(flow), "/")
library(tidyr); library(dplyr)
df_rel <- as.data.frame(flow_rel) %>% tibble::rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "group", values_to = "rel")
df_rel$group <- factor(df_rel$group, levels = names(color_rank))

ord <- rownames(flow)[order(rowSums(flow), decreasing = TRUE)]
df_rel$pathway <- factor(df_rel$pathway, levels = rev(ord))

ggplot(df_rel, aes(x = rel, y = pathway, fill = group)) +
  geom_col(width = 0.7, colour = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = 1/7, linetype = "dashed", colour = "grey35") +
  scale_fill_manual(values = color_rank) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 7.5),
        axis.text.x = element_text(colour = "black"),
        legend.title = element_blank()) +
  labs(x = "Relative information flow", y = NULL)

#####(3) focused on immune pathways
pw <- Reduce(union, lapply(names(chat_merge@netP), function(g) chat_merge@netP[[g]]$pathways))
length(pw)

flow <- sapply(names(chat_merge@netP), function(g) {
  x <- chat_merge@netP[[g]]
  v <- apply(x$prob, 3, sum); names(v) <- x$pathways
  v[pw]
})
flow[is.na(flow)] <- 0
rownames(flow) <- pw
det <- sapply(names(chat_merge@netP), function(g) pw %in% chat_merge@netP[[g]]$pathways)
rownames(det) <- pw
tab <- data.frame(pathway = pw, n_detected = rowSums(det), round(flow, 4), check.names = FALSE)
tab[order(tab$n_detected, decreasing = TRUE), ]


sel <- c("MHC-I","MHC-II","ICAM","ITGB2","SELPLG","CD99","LCK","CCL","CLEC","GALECTIN","IFN-II","ADGRE5","TNF","IL1","CD6","CD40","CD86","IL6","TGFb")
sel <- sel[sel %in% rownames(flow)]

flow_sel <- flow[sel, ]
flow_rel <- sweep(flow_sel, 1, rowSums(flow_sel), "/")

library(tidyr); library(dplyr)
df <- as.data.frame(flow_rel) %>% tibble::rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "group", values_to = "rel")
df$group   <- factor(df$group, levels = rev(names(color_rank)))
df$pathway <- factor(df$pathway, levels = rev(sel))

ggplot(df, aes(x = rel, y = pathway, fill = group)) +
  geom_col(width = 0.7, colour = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = 1/7, linetype = "dashed", colour = "grey35") +
  scale_fill_manual(values = color_rank) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 9),
        axis.text.x = element_text(colour = "black"),
        legend.title = element_blank()) +
  labs(x = "Relative information flow", y = NULL)

######7 differentially network (each group: polygenic and 5 gene mutations VS HC)
celltype_levels <- c("Naive CD4 T","Memory CD4 T",
                     "Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                     "γδT","Naive B","Memory B","Plasma","CD14 Monocyte",
                     "CD16 Monocyte","pDC","Platelet","Erythrocyte")
cell_color_vec <- c("#847AB3","#96AF95","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                    "#F1BAAB","#BB0021", "#D9C2D9","#FEC260","#AF478A",
                    "#EC706E","#B5BBE3","#E2A9C9","#CBDAA9")

plist_diff <- lapply(seq_along(chat_list)[-1], function(i) {
  nm <- paste0(names(chat_list)[i], " vs HC")
  as.ggplot(function() netVisual_diffInteraction(
    chat_merge, comparison = c(1, i), measure = "weight", weight.scale = TRUE,
    color.use = color_ct2, title.name = nm,
    vertex.weight = 6, vertex.size.max = 5, vertex.label.cex = 0.75,
    edge.width.max = 5, edge.curved = 0.2))
})
names(plist_diff) <- paste0(names(chat_list)[-1], "_vs_HC")
wrap_plots(plist_diff, ncol = 3)


######8 circors plot focused on immune pathway based on stacked result form above analysis for each group, seperately
sel <- c("MHC-I","MHC-II","ICAM","ITGB2","SELPLG","CD99","LCK","CCL","CLEC","GALECTIN","IFN-II","ADGRE5","TNF","IL1")
lab_full <- c(HC = "Healthy controls (C1-C6)", Poly = "Polygenic JIA (poly1-6)",
              S366Ter = "LACC1 (157LACC1)", I67R = "LPIN2 (182LPIN2)",
              M513K = "NOD2-M513K (181NOD2)", C495Y = "NOD2-C495Y (183NOD2)",
              E250K = "PSTPIP1-E250K (190+200)")

pdf("CC_chord_sel.pdf", width = 14, height = 10)
for (s in sel) {
  avail <- sapply(chat_list, function(x) s %in% x@netP$pathways)
  print(avail)
  grp_draw <- names(chat_list)[avail]
  if (!length(grp_draw)) next
  n_draw <- length(grp_draw); ncol_d <- min(4, n_draw); nrow_d <- ceiling(n_draw / ncol_d)
  circos.clear()
  par(mfrow = c(nrow_d, ncol_d), xpd = TRUE, mar = c(0.5, 0.5, 1.8, 0.5), cex = 0.45)
  for (g in grp_draw) {
    netVisual_aggregate(chat_list[[g]], signaling = s, layout = "chord", vertex.size = 6,
                        vertex.label.cex = 0.55, color.use = color_ct, title.name = lab_full[g])
  }
  par(mfrow = c(1, 1)); plot.new(); circos.clear()
}
par(mfrow = c(1, 1))
dev.off()


######9  patrial circos  
ct_keep <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
             "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC")

celltype_levels <- c("Naive CD4 T", "Memory CD4 T", "Naive CD8 T", "NK", "Treg",
                     "Cytotoxic CD8 T", "γδT", "Naive B", "Memory B", "Plasma",
                     "CD14 Monocyte", "CD16 Monocyte", "pDC", "Platelet", "Erythrocyte")
cell_color_vec <- c("#847AB3", "#96AF95", "#6699CC", "#86C7B4", "#007ABA", "#9CD2ED",
                    "#F1BAAB", "#BB0021", "#D9C2D9", "#FEC260", "#AF478A", "#EC706E",
                    "#B5BBE3", "#E2A9C9", "#CBDAA9")
color_ct2 <- setNames(cell_color_vec, celltype_levels)[ct_keep]
sel <- c("MHC-I","MHC-II","ICAM","ITGB2","SELPLG","CD99","LCK","CCL","CLEC","GALECTIN","IFN-II","ADGRE5","TNF","IL1","CD6","CD40","CD86","IL6","TGFb")

bipartite_net <- function(cc, path = NULL, ct = ct_keep, q = 0.90, ttl = NULL) {
  if (is.null(path)) {
    m <- cc@net$weight
  } else {
    idx <- which(cc@netP$pathways == path)
    if (!length(idx)) return(NULL)
    m <- cc@netP$prob[, , idx]
  }
  m <- m[intersect(rownames(m), ct), intersect(colnames(m), ct), drop = FALSE]
  thr <- quantile(m[m > 0], q, na.rm = TRUE)
  m[m < thr] <- 0
  ed <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(ed) <- c("s", "t", "w")
  ed <- ed[ed$w > 0 & ed$s != ed$t, ]
  if (!nrow(ed)) return(NULL)
  lib <- sort(unique(ed$s)); tib <- sort(unique(ed$t))
  nd <- data.frame(
    name  = c(paste0(lib, "__s"), paste0(tib, "__t")),
    label = c(lib, tib),
    side  = rep(c("s", "t"), c(length(lib), length(tib)))
  )
  yv <- setNames(seq(1, 0, length.out = length(union(lib, tib))), union(lib, tib))
  nd$y <- yv[nd$label]
  nd$x <- ifelse(nd$side == "s", 0, 1)
  ed$from <- paste0(ed$s, "__s"); ed$to <- paste0(ed$t, "__t")
  g <- graph_from_data_frame(ed[, c("from", "to", "w", "s")], vertices = nd, directed = TRUE)
  p <- ggraph(g, layout = "manual", x = nd$x, y = nd$y) +
    geom_edge_link(aes(edge_width = w, edge_colour = s),
                   alpha = 0.75, end_cap = circle(5, "mm"), start_cap = circle(5, "mm"),
                   arrow = arrow(angle = 18, length = unit(2, "mm"), type = "closed")) +
    geom_node_point(aes(filter = side == "s", colour = label), size = 4.5) +
    geom_node_point(aes(filter = side == "t", colour = label), size = 3.5, shape = 21, fill = NA, stroke = 1) +
    geom_node_text(aes(filter = side == "s", label = label), nudge_x = -0.05, hjust = 1, size = 2.5) +
    geom_node_text(aes(filter = side == "t", label = label), nudge_x = 0.05, hjust = 0, size = 2.5) +
    scale_edge_width(range = c(0.3, 2.6), guide = "none") +
    scale_edge_colour_manual(values = color_ct2, guide = "none") +
    scale_colour_manual(values = color_ct2, guide = "none") +
    coord_cartesian(xlim = c(-0.5, 1.5), clip = "off") +
    theme_void() + theme(plot.margin = margin(15, 20, 8, 20))
  if (!is.null(ttl)) p <- p + ggtitle(ttl)
  p
}

grp_all <- names(chat_list)
path_sets <- split(sel, ceiling(seq_along(sel) / 3))

pdf("CC_bipartite.pdf", width = 28, height = 13.5)
for (chunk in path_sets) {
  plots <- list()
  for (s in chunk) {
    ps <- lapply(grp_all, function(g) {
      t <- if (g == grp_all[1]) paste0(s, "  |  ", lab_full[g]) else lab_full[g]
      bipartite_net(chat_list[[g]], path = s, ttl = t)
    })
    ps <- ps[!sapply(ps, is.null)]
    plots <- c(plots, ps)
  }
  if (length(plots)) print(wrap_plots(plots, ncol = 7))
}
dev.off()


####check
edge_num <- function(g, path = "IL1", s = "pDC", t = c("CD14 Monocyte","CD16 Monocyte")){
  cc  <- chat_list[[g]]
  idx <- which(cc@netP$pathways == path)
  if (!length(idx)) return(NULL)
  m <- cc@netP$prob[, , idx]
  data.frame(group = g, target = t, prob = signif(m[s, t], 4),
             grp_max = signif(max(m), 3), grp_sum = signif(sum(m), 3), row.names = NULL)
}
do.call(rbind, lapply(names(chat_list), edge_num))

###change part 9 as the function
library(igraph)
library(ggraph)
library(patchwork)

##### CONFIG (already defined in your CC environment)
# ct_keep, color_ct2, lab_full, chat_list are expected to exist in the environment
ct_keep <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
             "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC")

celltype_levels <- c("Naive CD4 T", "Memory CD4 T", "Naive CD8 T", "NK", "Treg",
                     "Cytotoxic CD8 T", "γδT", "Naive B", "Memory B", "Plasma",
                     "CD14 Monocyte", "CD16 Monocyte", "pDC", "Platelet", "Erythrocyte")
cell_color_vec <- c("#847AB3", "#96AF95", "#6699CC", "#86C7B4", "#007ABA", "#9CD2ED",
                    "#F1BAAB", "#BB0021", "#D9C2D9", "#FEC260", "#AF478A", "#EC706E",
                    "#B5BBE3", "#E2A9C9", "#CBDAA9")
color_ct2 <- setNames(cell_color_vec, celltype_levels)[ct_keep]

options(future.globals.maxSize = 4 * 1024^3)
setwd("E:/Cohort PPT/JIA/code/CellCom")

short_lab <- c("HC", "Poly", "S366Ter", "I67R", "M513K", "C495Y", "E250K")
lab_full <- c(
  HC       = "Healthy controls (C1-C6)",
  Poly     = "Polygenic JIA (poly1-6)",
  S366Ter  = "LACC1 (157LACC1)",
  I67R     = "LPIN2 (182LPIN2)",
  M513K    = "NOD2-M513K (181NOD2)",
  C495Y    = "NOD2-C495Y (183NOD2)",
  E250K    = "PSTPIP1-E250K (190+200)"
)

chat_list <- list()
for (g in short_lab) {
  f <- paste0("chat_", g, ".RDS")
  if (file.exists(f)) {
    chat_list[[g]] <- readRDS(f)
    message("loaded: ", g)
  } else {
    message("missing: ", f)
  }
}

print(names(chat_list))
print(sapply(chat_list, nrow))
chat_merge<-readRDS("JIA_cellchat_mutation7_merge.RDS")

##### BIPARTITE NET: SINGLE PATHWAY x SINGLE CHAT OBJECT
bipartite_net <- function(cc, path = NULL, ct = ct_keep, q = 0.90, ttl = NULL) {
  if (is.null(path)) {
    m <- cc@net$weight
  } else {
    idx <- which(cc@netP$pathways == path)
    if (!length(idx)) return(NULL)
    m <- cc@netP$prob[, , idx]
  }
  m <- m[intersect(rownames(m), ct), intersect(colnames(m), ct), drop = FALSE]
  thr <- quantile(m[m > 0], q, na.rm = TRUE)
  m[m < thr] <- 0
  ed <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(ed) <- c("s", "t", "w")
  ed <- ed[ed$w > 0 & ed$s != ed$t, ]
  if (!nrow(ed)) return(NULL)
  lib <- sort(unique(ed$s)); tib <- sort(unique(ed$t))
  nd <- data.frame(
    name  = c(paste0(lib, "__s"), paste0(tib, "__t")),
    label = c(lib, tib),
    side  = rep(c("s", "t"), c(length(lib), length(tib)))
  )
  yv <- setNames(seq(1, 0, length.out = length(union(lib, tib))), union(lib, tib))
  nd$y <- yv[nd$label]
  nd$x <- ifelse(nd$side == "s", 0, 1)
  ed$from <- paste0(ed$s, "__s"); ed$to <- paste0(ed$t, "__t")
  g <- graph_from_data_frame(ed[, c("from", "to", "w", "s")], vertices = nd, directed = TRUE)
  p <- ggraph(g, layout = "manual", x = nd$x, y = nd$y) +
    geom_edge_link(aes(edge_width = w, edge_colour = s),
                   alpha = 0.75, end_cap = circle(5, "mm"), start_cap = circle(5, "mm"),
                   arrow = arrow(angle = 18, length = unit(2, "mm"), type = "closed")) +
    geom_node_point(aes(filter = side == "s", colour = label), size = 4.5) +
    geom_node_point(aes(filter = side == "t", colour = label), size = 3.5, shape = 21, fill = NA, stroke = 1) +
    geom_node_text(aes(filter = side == "s", label = label), nudge_x = -0.05, hjust = 1, size = 2.5) +
    geom_node_text(aes(filter = side == "t", label = label), nudge_x = 0.05, hjust = 0, size = 2.5) +
    scale_edge_width(range = c(0.3, 2.6), guide = "none") +
    scale_edge_colour_manual(values = color_ct2, guide = "none") +
    scale_colour_manual(values = color_ct2, guide = "none") +
    coord_cartesian(xlim = c(-0.5, 1.5), clip = "off") +
    theme_void() + theme(plot.margin = margin(15, 20, 8, 20))
  if (!is.null(ttl)) p <- p + ggtitle(ttl)
  p
}

bipartite_grid <- function(chat_list, pathways, groups = names(chat_list),
                           ct = ct_keep, q = 0.90, lab_full = NULL) {
  plots <- list()
  for (g in groups) {
    for (s in pathways) {
      ttl <- s
      if (!is.null(lab_full)) ttl <- paste0(lab_full[g], " | ", s)
      p <- bipartite_net(chat_list[[g]], path = s, ct = ct, q = q, ttl = ttl)
      if (!is.null(p)) plots[[length(plots) + 1]] <- p
    }
  }
  wrap_plots(plots, nrow = length(groups))
}

##### CALL EXAMPLE
pg <- bipartite_grid(
  chat_list = chat_list,
  pathways  = c("IFN-II", "IL1"),
  groups    = short_lab,
  lab_full  = lab_full
)
pg

####9.1 revised
library(igraph)
library(ggraph)
library(patchwork)

celltype_levels <- c("Naive CD4 T","Memory CD4 T","Naive CD8 T","NK","Treg","Cytotoxic CD8 T",
                     "γδT","Naive B","Memory B","Plasma","CD14 Monocyte","CD16 Monocyte","pDC",
                     "Platelet","Erythrocyte")
cell_color_vec <- c("#847AB3","#96AF95","#6699CC","#86C7B4","#007ABA","#9CD2ED",
                    "#F1BAAB","#BB0021","#D9C2D9","#FEC260","#AF478A","#EC706E",
                    "#B5BBE3","#E2A9C9","#CBDAA9")
color_ct2 <- setNames(cell_color_vec, celltype_levels)[ct_keep]

empty_panel <- function(msg, ttl = NULL){
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 3, colour = "grey45") +
    theme_void() + theme(plot.margin = margin(12, 16, 6, 16))
  if(!is.null(ttl)) p <- p + ggtitle(ttl) +
      theme(plot.title = element_text(size = 7, hjust = 0.5, lineheight = 0.95, colour = "black")) +
      theme(plot.margin = margin(4, 16, 6, 16))
  p
}

pw_stats <- function(path, ct = ct_keep){
  vals <- numeric(0); mx <- numeric(0)
  for(g in names(chat_list)){
    cc <- chat_list[[g]]; idx <- which(cc@netP$pathways == path)
    if(!length(idx)) next
    m <- cc@netP$prob[, , idx]
    m <- m[intersect(rownames(m), ct), intersect(colnames(m), ct), drop = FALSE]
    vals <- c(vals, m[m > 0]); mx <- c(mx, max(m, na.rm = TRUE))
  }
  list(vals = vals, max = if(length(mx)) max(mx, na.rm = TRUE) else 0)
}

bipartite_net <- function(cc, path = NULL, ct = ct_keep, thr = NULL, w_max = NULL,
                          ttl = NULL, label_size = 2.2){
  if(is.null(path)) m <- cc@net$weight else {
    idx <- which(cc@netP$pathways == path)
    if(!length(idx)) return(empty_panel("pathway absent", ttl))
    m <- cc@netP$prob[, , idx]
  }
  m <- m[intersect(rownames(m), ct), intersect(colnames(m), ct), drop = FALSE]
  tot <- sum(m, na.rm = TRUE)
  if(is.null(thr)) thr <- quantile(m[m > 0], 0.90, na.rm = TRUE)
  m[m < thr] <- 0
  ed <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(ed) <- c("s","t","w")
  ed <- ed[ed$w > 0 & ed$s != ed$t, ]
  if(!nrow(ed))
    return(empty_panel(paste0("all edges < ", signif(thr, 2), "\npathway total = ", signif(tot, 2)), ttl))
  if(is.null(w_max)) w_max <- max(ed$w)
  
  yv <- setNames(seq(1, 0, length.out = length(ct)), ct)
  ls <- unique(ed$s); lt <- unique(ed$t)
  nd <- rbind(
    data.frame(name = paste0(ls, "__s"), label = ls, side = "s", x = 0, y = unname(yv[ls]),
               stringsAsFactors = FALSE),
    data.frame(name = paste0(lt, "__t"), label = lt, side = "t", x = 1, y = unname(yv[lt]),
               stringsAsFactors = FALSE))
  ed$from <- paste0(ed$s, "__s"); ed$to <- paste0(ed$t, "__t")
  g <- graph_from_data_frame(ed[, c("from","to","w","s")], vertices = nd, directed = TRUE)
  
  p <- ggraph(g, layout = "manual", x = nd$x, y = nd$y) +
    geom_edge_link(aes(edge_width = w, edge_colour = s), alpha = 0.75,
                   end_cap = circle(5, "mm"), start_cap = circle(5, "mm"),
                   arrow = arrow(angle = 18, length = unit(2, "mm"), type = "closed")) +
    geom_node_point(aes(filter = side == "s", colour = label), size = 3.6) +
    geom_node_point(aes(filter = side == "t", colour = label), size = 3.0,
                    shape = 21, fill = NA, stroke = 0.8) +
    geom_node_text(aes(filter = side == "s", label = label), nudge_x = -0.05, hjust = 1,
                   size = label_size, colour = "grey20") +
    geom_node_text(aes(filter = side == "t", label = label), nudge_x = 0.05, hjust = 0,
                   size = label_size, colour = "grey20") +
    scale_edge_width(range = c(0.25, 2.4), limits = c(0, w_max), guide = "none") +
    scale_edge_colour_manual(values = color_ct2, guide = "none") +
    scale_colour_manual(values = color_ct2, guide = "none") +
    coord_cartesian(xlim = c(-0.42, 1.42), clip = "off") +
    theme_void() + theme(plot.margin = margin(4, 16, 6, 16))
  if(!is.null(ttl)) p <- p + ggtitle(ttl) +
    theme(plot.title = element_text(size = 7, hjust = 0.5, lineheight = 0.95, colour = "black"))
  p
}

sel_ok <- c("MHC-I","MHC-II","CLEC","CCL","GALECTIN","ICAM","ITGB2","SELPLG","ADGRE5","TGFb","CD6","IFN-II","IL6","IL1","TNF")
sel_ok <- sel_ok[sel_ok %in% unique(unlist(lapply(chat_list, function(cc) cc@netP$pathways)))]
floor_abs <- 1e-3

grp_all <- names(chat_list)
path_sets <- split(sel_ok, ceiling(seq_along(sel_ok) / 3))

pdf("CC_bipartite_globalThr.pdf", width = 30, height = 14)
for(chunk in path_sets){
  plots <- list()
  for(s in chunk){
    st <- pw_stats(s)
    thr_s <- max(quantile(st$vals, 0.90, na.rm = TRUE), floor_abs)
    plots <- c(plots, lapply(grp_all, function(g)
      bipartite_net(chat_list[[g]], path = s, thr = thr_s, w_max = st$max,
                    ttl = paste0(s, "\n", lab_full[g]))))
  }
  print(wrap_plots(plots, ncol = 7))
}
dev.off()


####check data
st    <- pw_stats("IL1")     
q90   <- quantile(st$vals, 0.90, na.rm = TRUE)   
thr_s <- max(q90, 1e-3)                         


#######9.2 circos changed
pw_max <- function(path, ct = ct_keep){
  mx <- 0
  for(g in names(chat_list)){
    cc <- chat_list[[g]]; idx <- which(cc@netP$pathways == path)
    if(!length(idx)) next
    m <- cc@netP$prob[, , idx]
    m <- m[intersect(rownames(m), ct), intersect(colnames(m), ct), drop = FALSE]
    mx <- max(mx, max(m, na.rm = TRUE))
  }
  mx
}

bipartite_net <- function(cc, path = NULL, ct = ct_keep, q = 0.90, w_max = NULL, ttl = NULL) {
  if (is.null(path)) {
    m <- cc@net$weight
  } else {
    idx <- which(cc@netP$pathways == path)
    if (!length(idx)) return(NULL)
    m <- cc@netP$prob[, , idx]
  }
  m <- m[intersect(rownames(m), ct), intersect(colnames(m), ct), drop = FALSE]
  thr <- quantile(m[m > 0], q, na.rm = TRUE)
  m[m < thr] <- 0
  ed <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  names(ed) <- c("s", "t", "w")
  ed <- ed[ed$w > 0 & ed$s != ed$t, ]
  if (!nrow(ed)) return(NULL)
  if (is.null(w_max) || !is.finite(w_max) || w_max <= 0) w_max <- max(ed$w)
  lib <- sort(unique(ed$s)); tib <- sort(unique(ed$t))
  nd <- data.frame(
    name  = c(paste0(lib, "__s"), paste0(tib, "__t")),
    label = c(lib, tib),
    side  = rep(c("s", "t"), c(length(lib), length(tib)))
  )
  yv <- setNames(seq(1, 0, length.out = length(union(lib, tib))), union(lib, tib))
  nd$y <- yv[nd$label]
  nd$x <- ifelse(nd$side == "s", 0, 1)
  ed$from <- paste0(ed$s, "__s"); ed$to <- paste0(ed$t, "__t")
  g <- graph_from_data_frame(ed[, c("from", "to", "w", "s")], vertices = nd, directed = TRUE)
  p <- ggraph(g, layout = "manual", x = nd$x, y = nd$y) +
    geom_edge_link(aes(edge_width = w, edge_colour = s),
                   alpha = 0.75, end_cap = circle(5, "mm"), start_cap = circle(5, "mm"),
                   arrow = arrow(angle = 18, length = unit(2, "mm"), type = "closed")) +
    geom_node_point(aes(filter = side == "s", colour = label), size = 4.5) +
    geom_node_point(aes(filter = side == "t", colour = label), size = 3.5, shape = 21, fill = NA, stroke = 1) +
    geom_node_text(aes(filter = side == "s", label = label), nudge_x = -0.05, hjust = 1, size = 2.5) +
    geom_node_text(aes(filter = side == "t", label = label), nudge_x = 0.05, hjust = 0, size = 2.5) +
    scale_edge_width(range = c(0.3, 2.6), limits = c(0, w_max), guide = "none") +
    scale_edge_colour_manual(values = color_ct2, guide = "none") +
    scale_colour_manual(values = color_ct2, guide = "none") +
    coord_cartesian(xlim = c(-0.5, 1.5), clip = "off") +
    theme_void() + theme(plot.margin = margin(15, 20, 8, 20))
  if (!is.null(ttl)) p <- p + ggtitle(ttl)
  p
}

grp_all <- names(chat_list)
path_sets <- split(sel, ceiling(seq_along(sel) / 3))

pdf("CC_bipartite.pdf", width = 28, height = 13.5)
for (chunk in path_sets) {
  plots <- list()
  for (s in chunk) {
    wm <- pw_max(s)
    ps <- lapply(grp_all, function(g) {
      bipartite_net(chat_list[[g]], path = s, w_max = wm,
                    ttl = paste0(s, "\n", lab_full[g]))
    })
    ps <- ps[!sapply(ps, is.null)]
    plots <- c(plots, ps)
  }
  if (length(plots)) print(wrap_plots(plots, ncol = 7))
}
dev.off()

