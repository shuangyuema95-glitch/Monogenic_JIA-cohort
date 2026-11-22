# -*- coding: utf-8 -*-
# CellPhoneDB v5.0.1 statistical analysis: permutation-based significance test
# run per genotype group (HC + Poly + 5 single-mutation groups)
# Input: per-group meta (cell ID - cell type) and counts (genes x cells, hgnc_symbol)
# Output: five result files per group (deconvoluted / deconvoluted_percents /
#         means / pvalues / significant_means)
# Run: conda activate cpdb && python 2_cpb.py
import os
import cellphonedb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

print(cellphonedb.version)

# Input directory holding the 14 txt files (7 meta + 7 counts)
DATA_DIR = "/public/home/yuxiaomingroup/msy/work/JIA"
# CellPhoneDB v5 background database zip (must point to the zip itself, not an unzipped folder)
DB_ZIP   = "/public/home/yuxiaomingroup/msy/work/cellphone/cell5.zip"
# Root output directory; one subdirectory per group is created below
OUT_ROOT = "/public/home/yuxiaomingroup/msy/work/cellphone/cpdb_out"

# Seven genotype groups, in cohort order
groups = ["HC", "Poly", "S366Ter", "I67R", "M513K", "C495Y", "E250K"]
# Prefixes of the five result files expected per group
five = ["deconvoluted", "deconvoluted_percents", "means", "pvalues", "significant_means"]

for g in groups:
    out = os.path.join(OUT_ROOT, g)
    os.makedirs(out, exist_ok=True)

    # Resume support: skip a group if all five result files already exist
    done = all(os.path.exists(os.path.join(out, f"statistical_analysis_{x}_{g}.txt")) for x in five)
    if done:
        print(f"=== {g} already done, skip ===", flush=True)
        continue

    # Per-group meta / counts absolute paths
    meta   = os.path.join(DATA_DIR, f"cpdb_{g}_meta.txt")
    counts = os.path.join(DATA_DIR, f"cpdb_{g}_counts.txt")
    print(f"=== running {g} ===", flush=True)
    print("meta exists:  ", os.path.exists(meta), meta, flush=True)
    print("counts exists:", os.path.exists(counts), counts, flush=True)

    # Statistical analysis: permutation test evaluates the significance of each ligand-receptor interaction
    cpdb_statistical_analysis_method.call(
        cpdb_file_path   = DB_ZIP,        # v5 background database zip
        meta_file_path   = meta,          # cell ID - cell type table
        counts_file_path = counts,        # genes x cells expression matrix (unnormalized counts)
        counts_data      = "hgnc_symbol", # gene identifier format: human HGNC symbol
        iterations       = 100,           # number of permutations; exploratory pass uses 100, raise to 1000 for final results
        threshold        = 0.1,           # minimum expression fraction of ligand/receptor within a cell type
        threads          = 8,             # parallel threads
        output_suffix    = g,             # suffix appended to output file names (group label)
        output_path      = out + "/"      # directory where result files are written
    )
    print("done:", g, flush=True)
