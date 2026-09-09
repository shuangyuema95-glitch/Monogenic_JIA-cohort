# ct2_run_all.py: CytoTRACE2 on all JIA pbmc1 cells (global, no subsetting)
import anndata
import pandas as pd
import os
import pickle
import gc
from cytotrace2_py.cytotrace2_py import cytotrace2

h5ad_path = "pbmc1.h5ad"
out_root  = "ct2_all"
os.makedirs(out_root, exist_ok=True)

print("Reading h5ad...")
adata = anndata.read_h5ad(h5ad_path)
raw_adata = adata.raw.to_adata()
raw_adata.var_names = adata.var_names.copy()
del adata
gc.collect()
print(f"  cells: {raw_adata.n_obs}, genes: {raw_adata.n_vars}")

print("CPM normalization...")
raw_counts = raw_adata.X.toarray()
lib_size   = raw_counts.sum(axis=1, keepdims=True)
cpm        = (raw_counts / lib_size) * 1e6
del raw_counts, lib_size
gc.collect()

expr_df = pd.DataFrame(cpm.T, index=raw_adata.var_names, columns=raw_adata.obs_names)
del cpm
gc.collect()

annot_df = pd.DataFrame({"phenotype": raw_adata.obs["celltype"].astype(str)},
                         index=raw_adata.obs_names)

expr_path = os.path.join(out_root, "expr_cpm.tsv")
annot_path = os.path.join(out_root, "annot.tsv")
print("Writing input TSVs...")
expr_df.to_csv(expr_path, sep="\t", index=True, header=True)
annot_df.to_csv(annot_path, sep="\t", index=True)
del expr_df, annot_df
gc.collect()

print("Running CytoTRACE2...")
res_df = cytotrace2(
    expr_path,
    annotation_path = annot_path,
    species         = "human",
    output_dir      = out_root,
    disable_plotting = True
)

for f in [expr_path, annot_path]:
    if os.path.exists(f):
        os.remove(f)

print("Joining results to obs...")
obs_full = raw_adata.obs.copy()
obs_full = obs_full.join(res_df)
out_pkl = os.path.join(out_root, "ct2_obs_all.pkl")
with open(out_pkl, "wb") as f:
    pickle.dump(obs_full, f)

print("=== DONE ===")
print(f"  pkl    : {out_pkl}")
print(f"  cells  : {obs_full.shape[0]}")
print(f"  columns: {obs_full.shape[1]}")

