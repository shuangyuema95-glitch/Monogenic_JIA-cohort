library(data.table)
library(ggplot2)

##### CONFIG
DATA_DIR <- "E:/Cohort PPT/JIA/code/CellCom/cpdb/cpdb_out"
OUT_DIR  <- "E:/Cohort PPT/JIA/code/CellCom"
groups   <- c("HC", "Poly", "S366Ter", "I67R", "M513K", "C495Y", "E250K")
P_CUT    <- 0.05

##### READ MEANS AND PVALUES PER GROUP
read_group <- function(g) {
  m <- fread(file.path(DATA_DIR, g, sprintf("statistical_analysis_means_%s.txt", g)), check.names = FALSE)
  p <- fread(file.path(DATA_DIR, g, sprintf("statistical_analysis_pvalues_%s.txt", g)), check.names = FALSE)
  m[, group := g]; p[, group := g]
  list(mean = m, pval = p)
}

to_long <- function(d, value_col) {
  pair_cols <- grep("\\|", names(d), value = TRUE)
  anno_cols <- setdiff(names(d), c(pair_cols, "group"))
  long <- melt(d, id.vars = c(anno_cols, "group"),
               measure.vars = pair_cols, variable.name = "pair", value.name = value_col)
  long[, c("cell1", "cell2") := tstrsplit(pair, "\\|")]
  long
}

##### BUILD FULL LONG TABLE (ALL LR PAIRS, SIGNIFICANT AND NOT)
res_list <- list()
for (g in groups) {
  rp <- read_group(g)
  lm <- to_long(rp$mean, "mean")
  lp <- to_long(rp$pval, "pval")
  anno_cols <- intersect(names(lm), names(lp))
  anno_cols <- setdiff(anno_cols, c("group", "pair", "cell1", "cell2", "mean", "pval"))
  key_cols  <- c(anno_cols, "group", "pair", "cell1", "cell2")
  dg <- merge(lm, lp[, c(key_cols, "pval"), with = FALSE], by = key_cols, all.x = TRUE)
  dg[, sig := !is.na(pval) & pval < P_CUT & !is.na(mean) & mean > 0]
  res_list[[g]] <- dg
}
res <- rbindlist(res_list, fill = TRUE)
fwrite(res, file.path(OUT_DIR, "cpdb_all_LR_long.csv"))
cat("all rows:", nrow(res), "| significant:", sum(res$sig, na.rm = TRUE), "\n")

##### SIGNIFICANT SUBSET TABLE
sig_dt <- res[sig == TRUE][order(group, cell1, cell2, pval)]
fwrite(sig_dt, file.path(OUT_DIR, "cpdb_significant_LR_long.csv"))

##### BUBBLE PLOT FUNCTION
plot_cpdb_bubble <- function(sig_dt, sender, receiver,
                             colors = c("#DCE9F5", "#6699CC"),
                             topn = 30, p_cut = P_CUT) {
  sub <- sig_dt[(cell1 == sender & cell2 == receiver) |
                  (cell1 == receiver & cell2 == sender)]
  
  agg <- sub[, .(mean = max(mean, na.rm = TRUE),
                 pval = min(pval, na.rm = TRUE)),
             by = .(group, interacting_pair)]
  
  top <- agg[, .(m = max(mean, na.rm = TRUE)), by = interacting_pair][
    order(-m)][1:min(topn, .N), interacting_pair]
  agg <- agg[interacting_pair %in% top]
  
  agg[, group := factor(group, levels = groups)]
  
  ggplot(agg, aes(x = group, y = reorder(interacting_pair, mean),
                  size = mean, color = mean)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(1, 5)) +
    scale_color_gradient(low = colors[1], high = colors[2], name = "mean") +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL,
         title = paste0(sender, "  <->  ", receiver),
         subtitle = paste0("significant: p < ", p_cut, " & mean > 0; ",
                           "top ", topn, " LR by max mean")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 8))
}

##### CALL EXAMPLE
p <- plot_cpdb_bubble(
  sig_dt,
  sender  = "CD16 Monocyte",
  receiver = "CD14 Monocyte",
  colors  = c("#DCE9F5", "#6699CC"),
  topn    = 30
)
p
