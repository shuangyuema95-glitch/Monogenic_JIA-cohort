###3D bar plots for GSEA results
########################################
#####multiple pathways merge into a figure
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

def plot_gsea_3d(filepath, box_size=0.3,
                 color_low="#B3D1E7", color_high="#F49D5C",
                 view=(22, -55), box_z=2.5):
    df = pd.read_csv(filepath, sep="\t")
    df["z"] = -np.log10(df["p.adjust"])

    agg = (df.sort_values("p.adjust")
             .groupby(["label", "mutant"], as_index=False)
             .first())

    labels = list(dict.fromkeys(agg["label"]))
    mutants = ["182LPIN2", "181NOD2", "183NOD2", "200PSTPIP1", "190PSTPIP1"]
    xi = {c: i for i, c in enumerate(labels)}
    yi = {m: i for i, m in enumerate(mutants)}

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection="3d")
    dx = dy = box_size

    cmap = LinearSegmentedColormap.from_list("nes_cmap", [color_low, color_high])
    norm = Normalize(vmin=agg["NES"].min(), vmax=agg["NES"].max())

    for _, r in agg.iterrows():
        x, y = xi[r["label"]], yi[r["mutant"]]
        ax.bar3d(x, y, 0, dx, dy, r["z"],
                 color=cmap(norm(r["NES"])),
                 edgecolor="k", linewidth=0.6, alpha=0.92)

    ax.set_xticks([i + dx/2 for i in range(len(labels))])
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_yticks([i + dy/2 for i in range(len(mutants))])
    ax.set_yticklabels(mutants)
    ax.set_zlabel("-log10(p.adjust)", labelpad=8)
    ax.set_box_aspect((len(labels), len(mutants), box_z))
    ax.view_init(elev=view[0], azim=view[1])

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label("NES")

    return fig, ax, agg


fig, ax, agg = plot_gsea_3d("E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\GSEA.txt",
                            color_low="#C8BBD5", color_high="#D6604D",box_size=0.38)
plt.show()



########################################
#### merge cell types, bar transparent
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

def plot_gsea_3d(filepath, box_size=0.4,
                 color_low="#B3D1E7", color_high="#F49D5C",
                 view=(22, -55), box_z=2.5, n_layers=100):
    rows = []
    with open(filepath) as f:
        f.readline()
        for line in f:
            t = line.split()
            if len(t) < 14:
                continue
            nes = float(t[4])
            padj = float(t[6])
            comparison = t[-1]
            mutant = t[-2]
            label = t[-3]
            j = len(t) - 4
            while j >= 0 and "/" not in t[j]:
                j -= 1
            celltype = " ".join(t[j+1:-3])
            rows.append((label, celltype, mutant, padj, nes))

    df = pd.DataFrame(rows, columns=["label", "celltype", "mutant", "padj", "NES"])
    agg = (df.sort_values("padj")
             .groupby(["label", "mutant"], as_index=False)
             .first())
    agg["z"] = -np.log10(agg["padj"])

    labels = list(dict.fromkeys(df["label"]))
    mutants = ["182LPIN2", "181NOD2", "183NOD2", "200PSTPIP1", "190PSTPIP1"]
    xi = {c: i for i, c in enumerate(labels)}
    yi = {m: i for i, m in enumerate(mutants)}

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection="3d")
    dx = dy = box_size

    cmap = LinearSegmentedColormap.from_list("nes_cmap", [color_low, color_high])
    norm = Normalize(vmin=agg["NES"].min(), vmax=agg["NES"].max())
    white = np.array([1.0, 1.0, 1.0, 1.0])

    for _, r in agg.iterrows():
        x, y = xi[r["label"]], yi[r["mutant"]]
        target = np.array(cmap(norm(r["NES"])))
        dz = r["z"] / n_layers
        for i in range(n_layers):
            frac = (i + 1) / n_layers
            color = (1 - frac) * white + frac * target
            ax.bar3d(x, y, i * dz, dx, dy, dz,
                     color=color, edgecolor="none")
        ax.bar3d(x, y, 0, dx, dy, r["z"],
                 facecolor=(0, 0, 0, 0), edgecolor="k", linewidth=0.6)

    ax.set_xticks([i + dx/2 for i in range(len(labels))])
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_yticks([i + dy/2 for i in range(len(mutants))])
    ax.set_yticklabels(mutants)
    ax.set_zlabel("-log10(p.adjust)", labelpad=8)
    ax.set_box_aspect((len(labels), len(mutants), box_z))
    ax.view_init(elev=view[0], azim=view[1])

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label("NES")

    return fig, ax, agg




fig, ax, agg = plot_gsea_3d("E:\\Cohort PPT\\JIA\\code\\CellPreprocess\\GSEA.txt",
                            color_low="#C8BBD5",color_high="#F49D5C")


########################################
####one pathway  one figure
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

PATHWAY_NAME = {
    "J": "IL6/JAK/STAT3",
    "A": "Interferon-alpha",
    "G": "Interferon-gamma",
    "N": "TNFA/NFKB",
}
PATHWAY_ORDER = ["J", "A", "G", "N"]

def plot_gsea_3d_facet(filepath, box_size=0.5,
                       color_low="#B3D1E7", color_high="#F49D5C",
                       view=(22, -55)):
    rows = []
    with open(filepath) as f:
        f.readline()
        for line in f:
            t = line.split()
            if len(t) < 14:
                continue
            nes = float(t[4])
            padj = float(t[6])
            mutant = t[-2]
            label = t[-3]
            j = len(t) - 4
            while j >= 0 and "/" not in t[j]:
                j -= 1
            celltype = " ".join(t[j+1:-3])
            rows.append((label, celltype, mutant, padj, nes))

    df = pd.DataFrame(rows, columns=["label", "celltype", "mutant", "padj", "NES"])
    df["suffix"] = df["label"].str.split("_").str[-1]
    df["ct"] = df["label"].str.split("_").str[0]
    df["z"] = -np.log10(df["padj"])

    mutants = list(dict.fromkeys(df["mutant"]))
    cts = list(dict.fromkeys(df["ct"]))
    yi = {m: i for i, m in enumerate(mutants)}
    xi = {c: i for i, c in enumerate(cts)}

    cmap = LinearSegmentedColormap.from_list("nes", [color_low, color_high])
    norm = Normalize(vmin=df["NES"].min(), vmax=df["NES"].max())

    fig = plt.figure(figsize=(18, 6))
    for k, suf in enumerate(PATHWAY_ORDER):
        sub = (df[df["suffix"] == suf]
               .sort_values("padj")
               .groupby(["ct", "mutant"], as_index=False)
               .first())
        ax = fig.add_subplot(1, 4, k+1, projection="3d")
        dx = dy = box_size
        for _, r in sub.iterrows():
            x = xi[r["ct"]]
            y = yi[r["mutant"]]
            ax.bar3d(x, y, 0, dx, dy, r["z"],
                     color=cmap(norm(r["NES"])),
                     edgecolor="k", linewidth=0.5, alpha=0.92)
            ax.text(x + dx/2, y + dy/2, r["z"] + 0.05,
                    f"{r['NES']:.2f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks([i + dx/2 for i in range(len(cts))])
        ax.set_xticklabels(cts, rotation=0)
        ax.set_yticks([i + dy/2 for i in range(len(mutants))])
        if k > 0:
            ax.set_yticklabels([])
        ax.set_zlabel("-log10(p.adjust)" if k == 0 else "")
        ax.set_title(PATHWAY_NAME[suf], fontsize=12)
        ax.set_box_aspect((len(cts), len(mutants), 2.5))
        ax.view_init(elev=view[0], azim=view[1])

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=fig.axes, shrink=0.5, pad=0.02)
    cbar.set_label("NES")

    return fig, df


fig, df = plot_gsea_3d_facet("C:/Users/Administrator/Desktop/GSEA.txt")
plt.show()


#####################################
####################################system prop 3D bar
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_clinical_3d(filepath, system_order, pathway_order,
                     pathway_colors, box_size=0.6, gap=0.2,
                     view=(25, -55), box_z=3):
    df = pd.read_csv(filepath, sep="\t")
    xi = {s: i for i, s in enumerate(system_order)}
    yi = {p: i for i, p in enumerate(pathway_order)}
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111, projection="3d")
    dx = dy = box_size
    offset = box_size / 2 + gap
    for _, r in df.iterrows():
        x = xi[r["System"]] * (box_size + gap)
        y = yi[r["Pathway"]] * (box_size + gap)
        ax.bar3d(x, y, 0, dx, dy, r["Percentage"],
                 color=pathway_colors[r["Pathway"]],
                 edgecolor="k", linewidth=0.5, alpha=1.0)
    ax.set_xticks([i * (box_size + gap) + box_size/2 for i in range(len(system_order))])
    ax.set_xticklabels(system_order, rotation=35, ha="right", fontsize=7)
    ax.set_yticks([i * (box_size + gap) + box_size/2 for i in range(len(pathway_order))])
    ax.set_yticklabels(pathway_order, fontsize=8)
    ax.set_zlabel("Percentage (%)", labelpad=8)
    ax.set_box_aspect((len(system_order)*(box_size+gap),
                       len(pathway_order)*(box_size+gap), box_z))
    ax.view_init(elev=view[0], azim=view[1])
    plt.tight_layout()
    return fig, ax


system_order = [
    "Growth and musculoskeletal system",
    "Immune system",
    "Systemic Manifestations",
    "Cutaneous-Mucosal system",
    "Hematologic system",
    "Digestive system",
    "Hepatosplenic system",
    "Endocrine system",
    "Respiratory system",
    "Nervous system",
    "Urinary system",
    "Cardiovascular system"
]

pathway_order = ["Interferon", "NFKB", "Immune metabolism", "Inflammasome", "Uncategoried", "Cell death"]
pathway_colors = {"NFKB": "#B3D1E7",
                  "Inflammasome": "#6699CC",
                  "Immune metabolism": "#847AB3",
                  "Interferon": "#EC706E",
                  "Cell death": "#FEC260",
                  "Uncategoried": "#EBB1A4"}

fig, ax = plot_clinical_3d(
    filepath="E:/Cohort PPT/JIA/code/CellPreprocess/clinical_3d.txt",
    system_order=system_order,
    pathway_order=pathway_order,
    pathway_colors=pathway_colors,
    box_size=0.25,
    gap=0.6,
    view=(30, -70),
    box_z=2.5
)
plt.show()



################加上渐变 不描边
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from mpl_toolkits.mplot3d import Axes3D

def gradient_bar3d(ax, x, y, z_bottom, dx, dy, height, cmap):
    n_steps = 20
    dz = height / n_steps
    for i in range(n_steps):
        z = z_bottom + i * dz
        frac = i / n_steps
        color = cmap(frac)
        ax.bar3d(x, y, z, dx, dy, dz, color=color, edgecolor="none", shade=False)

def plot_clinical_3d(filepath, system_order, pathway_order,
                     pathway_colors, box_size=0.25, gap=0.6,
                     view=(30, -70), box_z=2):
    df = pd.read_csv(filepath, sep="\t")
    xi = {s: i for i, s in enumerate(system_order)}
    yi = {p: i for i, p in enumerate(pathway_order)}
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111, projection="3d")
    dx = dy = box_size
    for _, r in df.iterrows():
        x = xi[r["System"]] * (box_size + gap)
        y = yi[r["Pathway"]] * (box_size + gap)
        base = np.array(to_rgb(pathway_colors[r["Pathway"]]))
        light = tuple(base * 0.4 + 0.6)
        dark = tuple(base)
        cmap = LinearSegmentedColormap.from_list("grad", [light, dark])
        gradient_bar3d(ax, x, y, 0, dx, dy, r["Percentage"], cmap)
    ax.set_xticks([i * (box_size + gap) + box_size/2 for i in range(len(system_order))])
    ax.set_xticklabels(system_order, rotation=35, ha="right", fontsize=7)
    ax.set_yticks([i * (box_size + gap) + box_size/2 for i in range(len(pathway_order))])
    ax.set_yticklabels(pathway_order, fontsize=8)
    ax.set_zlabel("Percentage (%)", labelpad=8)
    ax.set_box_aspect((len(system_order)*(box_size+gap),
                       len(pathway_order)*(box_size+gap), box_z))
    ax.view_init(elev=view[0], azim=view[1])
    plt.tight_layout()
    return fig, ax

pathway_colors = {"NFKB": "#B3D1E7",
                  "Inflammasome": "#0B71AB",
                  "Immune metabolism": "#847AB3",
                  "Interferon": "#EC706E",
                  "Cell death": "#FEC260",
                  "Uncategoried": "#EBB1A4"}

system_order = [
    "Cutaneous-Mucosal system",
    "Growth and musculoskeletal system",
    "Immune system",
    "Systemic Manifestations",
    "Hematologic system",
    "Digestive system",
    "Respiratory system",
    "Endocrine system",
    "Hepatosplenic system",
    "Cardiovascular system",
    "Urinary system",
    "Nervous system"
]
pathway_order = ["Interferon", "NFKB", "Immune metabolism", "Inflammasome", "Uncategoried", "Cell death"]

fig, ax = plot_clinical_3d(
    filepath="E:/Cohort PPT/JIA/code/CellPreprocess/clinical_3d.txt",
    system_order=system_order,
    pathway_order=pathway_order,
    pathway_colors=pathway_colors,
    box_size=0.25,
    gap=0.6,
    view=(30, -70),
    box_z=2
)
plt.show()


#######加上渐变 描边
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from mpl_toolkits.mplot3d import Axes3D

def gradient_bar3d(ax, x, y, z_bottom, dx, dy, height, cmap):
    n_steps = 20
    dz = height / n_steps
    for i in range(n_steps):
        z = z_bottom + i * dz
        frac = i / n_steps
        color = cmap(frac)
        ax.bar3d(x, y, z, dx, dy, dz, color=color, edgecolor="none", shade=False)
    ax.bar3d(x, y, 0, dx, dy, height, color=(0,0,0,0), edgecolor="#888888", linewidth=0.3, shade=False)

def plot_clinical_3d(filepath, system_order, pathway_order,
                     pathway_colors, box_size=0.25, gap=0.6,
                     view=(30, -70), box_z=2):
    df = pd.read_csv(filepath, sep="\t")
    xi = {s: i for i, s in enumerate(system_order)}
    yi = {p: i for i, p in enumerate(pathway_order)}
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(111, projection="3d")
    dx = dy = box_size
    for _, r in df.iterrows():
        x = xi[r["System"]] * (box_size + gap)
        y = yi[r["Pathway"]] * (box_size + gap)
        base = np.array(to_rgb(pathway_colors[r["Pathway"]]))
        light = tuple(base * 0.4 + 0.6)
        dark = tuple(base)
        cmap = LinearSegmentedColormap.from_list("grad", [light, dark])
        gradient_bar3d(ax, x, y, 0, dx, dy, r["Percentage"], cmap)
    ax.set_xticks([i * (box_size + gap) + box_size/2 for i in range(len(system_order))])
    ax.set_xticklabels(system_order, rotation=35, ha="right", fontsize=7)
    ax.set_yticks([i * (box_size + gap) + box_size/2 for i in range(len(pathway_order))])
    ax.set_yticklabels(pathway_order, fontsize=8)
    ax.set_zlabel("Percentage (%)", labelpad=8)
    ax.set_box_aspect((len(system_order)*(box_size+gap),
                       len(pathway_order)*(box_size+gap), box_z))
    ax.view_init(elev=view[0], azim=view[1])
    plt.tight_layout()
    return fig, ax

pathway_colors = {"NFKB": "#B3D1E7",
                  "Inflammasome": "#0B71AB",
                  "Immune metabolism": "#847AB3",
                  "Interferon": "#EC706E",
                  "Cell death": "#FEC260",
                  "Uncategoried": "#EBB1A4"}

system_order = [
    "Growth and musculoskeletal system",
    "Immune system",
    "Systemic Manifestations",
    "Cutaneous-Mucosal system",
    "Hematologic system",
    "Digestive system",
    "Hepatosplenic system",
    "Endocrine system",
    "Respiratory system",
    "Nervous system",
    "Urinary system",
    "Cardiovascular system"
]
pathway_order = ["Interferon", "NFKB", "Immune metabolism", "Inflammasome", "Uncategoried", "Cell death"]

fig, ax = plot_clinical_3d(
    filepath="E:/Cohort PPT/JIA/code/CellPreprocess/clinical_3d.txt",
    system_order=system_order,
    pathway_order=pathway_order,
    pathway_colors=pathway_colors,
    box_size=0.25,
    gap=0.6,
    view=(30, -70),
    box_z=2.5
)
plt.show()
