import os
import numpy as np
import pandas as pd
import scanpy as sc
import palantir
import matplotlib.pyplot as plt
from scipy import sparse

os.chdir("E:\\Cohort PPT\\JIA\\code\\CellTrajectory")

# ---------- 1. READ + SUBSET ----------
adata = sc.read("pbmc1.h5ad")
print("All cells:", adata.n_obs)
subset = adata[adata.obs['celltype'].isin(['CD14 Monocyte', 'CD16 Monocyte'])].copy()
print("Monocytes:", subset.n_obs)
print(subset.obs['celltype'].value_counts())

# ---------- 2. PREPROCESS ----------
sc.pp.normalize_per_cell(subset)
subset.X = subset.X.toarray()
palantir.preprocess.log_transform(subset)
sc.pp.highly_variable_genes(subset, n_top_genes=1500, flavor='cell_ranger')
sc.pp.pca(subset)
pca_projections = pd.DataFrame(subset.obsm['X_pca'], index=subset.obs_names)
umap = pd.DataFrame(subset.obsm['X_umap'] if 'X_umap' in subset.obsm else subset.obsm['X_pca'][:, :2],
                    index=subset.obs_names)
umap.columns = ["x", "y"]

# ---------- 3. DIFFUSION MAPS + MAGIC ----------
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
subset.X = sparse.csr_matrix(subset.X)
subset.layers['MAGIC_imputed_data'] = palantir.utils.run_magic_imputation(subset, dm_res)

# ---------- 4. START CELL (highest CD14) ----------
cd14_expr = subset[:, 'CD14'].X
cd14_expr = cd14_expr.toarray().flatten() if hasattr(cd14_expr, 'toarray') else cd14_expr.flatten()
start_cell = subset.obs_names[np.argmax(cd14_expr)]
print("Start cell:", start_cell, "| celltype:", subset.obs.loc[start_cell, 'celltype'])

# ---------- 5. TERMINAL STATE (highest CD16) ----------
cd16_expr = subset[:, 'FCGR3A'].X
cd16_expr = cd16_expr.toarray().flatten() if hasattr(cd16_expr, 'toarray') else cd16_expr.flatten()
term_cell = subset.obs_names[np.argmax(cd16_expr)]
print("Terminal cell:", term_cell, "| celltype:", subset.obs.loc[term_cell, 'celltype'])

# ---------- 6. RUN PALANTIR ----------
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500,
                                    terminal_states=[term_cell])
pr_res.branch_probs.columns = ['CD16']



# ---------- 7. EXPORT ----------
pr_res.pseudotime.to_csv("palantir_pseudotime.csv")
