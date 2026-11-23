# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

CPDB_DIR = Path("/public/home/yuxiaomingroup/msy/work/cellphone/cpdb_out")
OUT_DIR  = Path("/public/home/yuxiaomingroup/msy/work/cellphone/cpdb_fig")
OUT_DIR.mkdir(exist_ok=True)

groups = ["HC", "Poly", "S366Ter", "I67R", "M513K", "C495Y", "E250K"]
P_CUT  = 0.05

ct_order = ["Naive CD4 T", "Memory CD4 T", "Naive CD8 T", "NK", "Treg",
            "Cytotoxic CD8 T", "γδT", "Naive B", "Memory B", "Plasma",
            "CD14 Monocyte", "CD16 Monocyte", "pDC"]

def count_matrix(g):
    p = pd.read_csv(CPDB_DIR / g / f"statistical_analysis_pvalues_{g}.txt", sep="\t")
    pair_cols = [c for c in p.columns if "|" in c]
    sig = (p[pair_cols] < P_CUT).sum(axis=0)
    mat = pd.DataFrame(0, index=ct_order, columns=ct_order, dtype=float)
    for col in pair_cols:
        a, b = col.split("|")
        if a in mat.index and b in mat.columns:
            mat.loc[a, b] = sig[col]
    return mat

mats = {g: count_matrix(g) for g in groups}
vmax = max(m.values.max() for m in mats.values())

fig, axes = plt.subplots(1, 7, figsize=(24, 6), sharey=True)
for i, (g, ax) in enumerate(zip(groups, axes)):
    sns.heatmap(
        mats[g], ax=ax, cmap="Reds", vmin=0, vmax=vmax,
        cbar=(i == 6), cbar_kws={"label": "sig. LR pairs"},
        xticklabels=True, yticklabels=(i == 0), square=True
    )
    ax.set_title(g, fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.tick_params(axis="y", labelsize=7)

fig.supylabel("Cell type", x=0.01)
plt.tight_layout()
plt.savefig(OUT_DIR / "cpdb_heatmap_7groups.pdf", bbox_inches="tight")
plt.close()
print("saved:", OUT_DIR / "cpdb_heatmap_7groups.pdf")

