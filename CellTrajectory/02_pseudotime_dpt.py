# ============================================
# JIA: PAGA + DPT on CD14/CD16 Monocytes
# ============================================
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# ---------- 0. READ + SUBSET ----------
adata = sc.read("pbmc1.h5ad")
adata.obs['celltype'] = adata.obs['celltype'].astype('category')
print("All cells:", adata.n_obs)
print(adata.obs['celltype'].value_counts())

subset = adata[adata.obs['celltype'].isin(['CD14 Monocyte', 'CD16 Monocyte'])].copy()
print("\nMonocytes:", subset.n_obs)
print(subset.obs['celltype'].value_counts())

# ---------- 0b. SAVE CD14 EXPRESSION BEFORE HVG FILTERING ----------
cd14_x = subset[:, 'CD14'].X
subset.obs['CD14_expr'] = cd14_x.toarray().flatten() if hasattr(cd14_x, 'toarray') else cd14_x.flatten()
print("CD14 expr range:", subset.obs['CD14_expr'].min(), "-", subset.obs['CD14_expr'].max())

# ---------- 1. PREPROCESS (within monocytes) ----------
sc.pp.highly_variable_genes(subset, n_top_genes=2000)
subset = subset[:, subset.var.highly_variable].copy()
sc.pp.scale(subset, max_value=10)
sc.tl.pca(subset, n_comps=30)
sc.pp.neighbors(subset, n_pcs=20, n_neighbors=15)
sc.tl.umap(subset)

# ---------- 2. DIFFUSION MAP ----------
sc.tl.diffmap(subset, n_comps=20)
print("\nDiffmap shape:", subset.obsm['X_diffmap'].shape)
sc.pl.diffmap(subset, color='celltype', use_raw=False, save='_monocyte.png')

# ---------- 3. PAGA (trajectory topology) ----------
sc.tl.paga(subset, groups='celltype')
sc.pl.paga(subset, save='_monocyte.png')

# ---------- 4. SET ROOT (CD14-highest cell, from obs) ----------
root_idx = int(subset.obs['CD14_expr'].argmax())
root_cell = subset.obs_names[root_idx]
print(f"\nRoot cell (highest CD14): {root_cell}, index {root_idx}")
print(f"Root celltype: {subset.obs.loc[root_cell, 'celltype']}")
subset.uns['iroot'] = root_idx

# ---------- 5. DPT ----------
sc.tl.dpt(subset)
print("\nDPT range:", subset.obs['dpt_pseudotime'].min(), "-", subset.obs['dpt_pseudotime'].max())
print(subset.obs[['celltype', 'dpt_pseudotime']].groupby('celltype').describe())

sc.pl.umap(subset, color=['dpt_pseudotime', 'celltype'], save='_dpt.png')

# ---------- 6. DENSITY PLOT ----------
g = sns.FacetGrid(
    subset.obs, row="celltype", aspect=4, height=2,
    sharex=True, sharey=False, row_order=['CD14 Monocyte', 'CD16 Monocyte']
)
g.map(sns.kdeplot, "dpt_pseudotime", fill=True, color="skyblue", alpha=0.7)
g.set_titles("{row_name}")
g.set_axis_labels("DPT", "Density")
g.fig.subplots_adjust(hspace=0.5)
plt.savefig("dpt_density.png", dpi=150, bbox_inches='tight')
plt.show()

# ---------- 7. EXPORT ----------
subset.obs[['celltype', 'dpt_pseudotime']].to_csv("dpt_pseudotime.csv", index=True)
print("\n=== DONE ===")
print("Saved: dpt_pseudotime.csv")

