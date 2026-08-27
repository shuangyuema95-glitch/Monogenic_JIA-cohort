import warnings
warnings.filterwarnings("ignore")
import os
import pandas as pd
import matplotlib.pyplot as plt
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
import arviz as az
import numpy as np

input = "/public/home/yuxiaomingroup/msy/work/JIA"
os.chdir(input)

def get_min_variance_ref_cell(df_count):
    count_mat = df_count.copy()
    group_col = "sample"
    cell_cols = [col for col in count_mat.columns if col != group_col]
    var_dict = {}
    for cell in cell_cols:
        group_mean = count_mat.groupby(group_col)[cell].mean()
        var_dict[cell] = np.var(group_mean)
    var_df = pd.DataFrame(list(var_dict.items()), columns=["cell_type", "group_variance"])
    var_df = var_df.sort_values("group_variance", ascending=True).reset_index(drop=True)
    ref_cell = var_df.iloc[0]["cell_type"]
    print(f"Auto selected reference cell (min group variance): {ref_cell}")
    print("Variance ranking (top5 least variable):")
    print(var_df.head(5))
    return ref_cell


def run_scCODA_contrast(
        csv_path,
        ref_group_label,
        target_groups,
        fdr=0.4,
        sampler="nuts",
        output_csv=""
):
    cell_counts = pd.read_csv(csv_path)
    cell_counts.columns = cell_counts.columns.str.strip()
    print("Cleaned columns:", cell_counts.columns.tolist())

    ref_cell = get_min_variance_ref_cell(cell_counts)

    data_all = dat.from_pandas(cell_counts, covariate_columns=["sample"])

    hc_samples = {"C1", "C2", "C3", "C4", "C5", "C6"}
    poly_samples = {"poly1", "poly2", "poly3", "poly4", "poly5", "poly6"}
    mono_samples = {"182LPIN2", "181NOD2", "183NOD2", "200PSTPIP1", "190PSTPIP1"}

    def map_broad(s):
        if s in hc_samples:
            return "HC"
        elif s in poly_samples:
            return "polygenic"
        elif s in mono_samples:
            return "monogenic"
        else:
            return np.nan

    data_all.obs["broad_group"] = data_all.obs["sample"].apply(map_broad)
    data_all = data_all[~data_all.obs["broad_group"].isna()].copy()

    all_res = []
    for group in target_groups:
        print(f"\n=== Running contrast: {group} VS {ref_group_label} | Sampler: {sampler} ===")
        sub_data = data_all[data_all.obs["broad_group"].isin([ref_group_label, group])]

        formula = f"C(broad_group, Treatment(reference='{ref_group_label}'))"
        print(f"Model formula: {formula}")

        model = mod.CompositionalAnalysis(
            sub_data,
            formula=formula,
            reference_cell_type=ref_cell
        )

        if sampler.lower() == "nuts":
            sim = model.sample_nuts()
        else:
            sim = model.sample_hmc()

        sim.set_fdr(est_fdr=fdr)
        eff = sim.effect_df.reset_index()
        cred = sim.credible_effects().reset_index(name="Significant")
        df = pd.merge(eff, cred, on=["Covariate", "Cell Type"])
        df["Contrast"] = f"{group}_vs_{ref_group_label}"
        all_res.append(df)

    final_df = pd.concat(all_res, ignore_index=True)
    final_df = final_df.rename(columns={"log2-fold change": "log2FC"})
    final_df.to_csv(output_csv, index=False)
    print(f"\nDone! Output file: {output_csv}")
    return final_df


# 1) polygenic VS HC
res1 = run_scCODA_contrast(
    csv_path="resultCount.csv",
    ref_group_label="HC",
    target_groups=["polygenic"],
    output_csv="scCODA_poly_vs_HC.csv"
)

# 2) monogenic VS HC
res2 = run_scCODA_contrast(
    csv_path="resultCount.csv",
    ref_group_label="HC",
    target_groups=["monogenic"],
    output_csv="scCODA_mono_vs_HC.csv"
)

# 3) monogenic VS polygenic
res3 = run_scCODA_contrast(
    csv_path="resultCount.csv",
    ref_group_label="polygenic",
    target_groups=["monogenic"],
    output_csv="scCODA_mono_vs_poly.csv"
)
