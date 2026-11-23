# -*- coding: utf-8 -*-
# CellPhoneDB dot/bubble plot: ligand-receptor pairs between two cell types
# point size = mean expression, color/significance from permutation p-values
import anndata as ad
import pandas as pd
import ktplotspy as kpy
from pathlib import Path

H5AD     = "/public/home/yuxiaomingroup/msy/work/JIA/pbmc1.h5ad"
CPDB_DIR = Path("/public/home/yuxiaomingroup/msy/work/cellphone/cpdb_out")
OUT_DIR  = Path("/public/home/yuxiaomingroup/msy/work/cellphone/cpdb_fig")
OUT_DIR.mkdir(exist_ok=True)

GROUP = "M513K"
CT1   = "CD16 Monocyte"
CT2   = "CD14 Monocyte"
# None = all significant pairs for this cell pair; or pass a gene list, e.g.
# GENES = ["CD44", "TYROBP", "ANXA1", "FPR1", "SELPLG", "ICAM1", "ITGB2"]
GENES = None

adata = ad.read_h5ad(H5AD)
print("obs columns:", adata.obs.columns.tolist())
print("celltypes:", adata.obs["celltype"].unique().tolist())

means = pd.read_csv(CPDB_DIR / GROUP / f"statistical_analysis_means_{GROUP}.txt", sep="\t")
pvals = pd.read_csv(CPDB_DIR / GROUP / f"statistical_analysis_pvalues_{GROUP}.txt", sep="\t")

p = kpy.plot_cpdb(
    adata=adata,
    cell_type1=CT1,
    cell_type2=CT2,
    means=means,
    pvals=pvals,
    celltype_key="celltype",
    genes=GENES,
    figsize=(13, 5),
    title=f"{GROUP}: {CT1} vs {CT2}",
    keep_id_cp_interaction=True
)

out_pdf = OUT_DIR / f"bubble_{GROUP}.pdf"
p.save(str(out_pdf))
print("saved:", out_pdf)

