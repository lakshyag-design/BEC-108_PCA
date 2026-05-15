"""
PCA on Breast Cancer Gene Expression Data (GSE5325)
Reproduces Figure 1a, 1b, and Figure 1c from the Nature Primer (Ringnér 2008).

Dependencies:
    pip install pandas numpy matplotlib scikit-learn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ─────────────────────────────────────────────
# 1.  LOAD DATA
# ─────────────────────────────────────────────

print("Loading data...")

# Gene expression matrix  (rows = patients, columns = gene IDs)
expr = pd.read_csv("data/filtered.tsv.gz", sep="\t", index_col=None, header=0)
# Strip whitespace from column names (original file has leading spaces)
expr.columns = expr.columns.str.strip()

# Class labels  (105 rows, no header; 1 = ER+, 0 = ER-)
labels = pd.read_csv("data/class.tsv", sep="\t", header=None, names=["label"])

# Gene ID → symbol mapping
gene_map = pd.read_csv(
    "data/columns.tsv.gz",
    sep="\t",
    comment="#",
    usecols=["ID", "GeneSymbol"],
)
gene_map["ID"] = gene_map["ID"].astype(str).str.strip()
# Keep first occurrence of each ID (some IDs appear twice with different annotations)
gene_map = gene_map.drop_duplicates(subset="ID", keep="first")
id_to_symbol = dict(zip(gene_map["ID"], gene_map["GeneSymbol"]))

print(f"Expression matrix shape: {expr.shape}  (patients × genes)")
print(f"Labels: {labels['label'].value_counts().to_dict()}")

# Confirm shapes match
assert expr.shape[0] == len(labels), "Row count mismatch between expr and labels!"

# ─────────────────────────────────────────────
# 2.  EXTRACT XBP1 AND GATA3
# ─────────────────────────────────────────────

XBP1_ID  = "4404"   # confirmed from columns.tsv.gz
GATA3_ID = "4359"

xbp1  = expr[XBP1_ID].values
gata3 = expr[GATA3_ID].values
label = labels["label"].values          # 1 = ER+, 0 = ER-

colors = np.where(label == 1, "red", "black")

# ─────────────────────────────────────────────
# 3.  FIGURE 1a  –  Scatter plot XBP1 vs GATA3
# ─────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(5, 5))

ax.scatter(gata3, xbp1, c=colors, s=25, marker="s", linewidths=0)

ax.set_xlabel("GATA3", fontstyle="italic", fontsize=13)
ax.set_ylabel("XBP1",  fontstyle="italic", fontsize=13)

er_pos = mpatches.Patch(color="red",   label="ER⁺")
er_neg = mpatches.Patch(color="black", label="ER⁻")
ax.legend(handles=[er_pos, er_neg], frameon=False, fontsize=11)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("figure1a_scatter.png", dpi=150)
plt.close()
print("Saved → figure1a_scatter.png")

# ─────────────────────────────────────────────
# 4.  FIGURE 1b  –  Scatter plot WITH PC1/PC2 arrows
#     PCA is done only on XBP1 & GATA3 (2-gene space)
# ─────────────────────────────────────────────

# PCA on just the two genes to get the principal axes in this 2D space
two_gene = np.column_stack([gata3, xbp1])   # shape (105, 2)
pca_2d = PCA(n_components=2)
pca_2d.fit(two_gene)

# Center of the data (origin of arrows)
cx, cy = np.mean(gata3), np.mean(xbp1)

# Scale arrows to be visible (proportional to explained variance)
scale = 2.0
pc1_vec = pca_2d.components_[0] * scale   # [gata3_component, xbp1_component]
pc2_vec = pca_2d.components_[1] * scale

fig, ax = plt.subplots(figsize=(5, 5))

ax.scatter(gata3, xbp1, c=colors, s=25, marker="s", linewidths=0)

# Draw PC1 arrow (extended both directions as a full line)
arrow_kw = dict(head_width=0.15, head_length=0.15,
                fc="black", ec="black", length_includes_head=True)

# PC1 arrow
ax.annotate("", xy=(cx + pc1_vec[0], cy + pc1_vec[1]),
            xytext=(cx - pc1_vec[0], cy - pc1_vec[1]),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5))
ax.text(cx + pc1_vec[0] + 0.15, cy + pc1_vec[1], "PC1", fontsize=11, fontweight="bold")

# PC2 arrow
ax.annotate("", xy=(cx + pc2_vec[0], cy + pc2_vec[1]),
            xytext=(cx - pc2_vec[0], cy - pc2_vec[1]),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5))
ax.text(cx + pc2_vec[0] + 0.1, cy + pc2_vec[1] + 0.1, "PC2", fontsize=11, fontweight="bold")

ax.set_xlabel("GATA3", fontstyle="italic", fontsize=13)
ax.set_ylabel("XBP1",  fontstyle="italic", fontsize=13)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("figure1b_scatter_pca_arrows.png", dpi=150)
plt.close()
print("Saved → figure1b_scatter_pca_arrows.png")

# ─────────────────────────────────────────────
# 6.  PCA on FULL expression matrix
# ─────────────────────────────────────────────

print("\nRunning PCA on full expression matrix...")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr.values)   # shape: (105, n_genes)

pca = PCA(n_components=min(100, expr.shape[0]))  # keep up to 100 PCs
X_pca = pca.fit_transform(X_scaled)             # shape: (105, n_components)

pc1_scores = X_pca[:, 0]   # projection of each patient onto PC1

print(f"Variance explained by PC1: {pca.explained_variance_ratio_[0]*100:.1f}%")
print(f"Variance explained by PC2: {pca.explained_variance_ratio_[1]*100:.1f}%")

# ─────────────────────────────────────────────
# 5.  FIGURE 1c  –  Projection onto PC1
#     Three rows: All / ER- / ER+
# ─────────────────────────────────────────────

er_pos_mask = label == 1
er_neg_mask = label == 0

fig, ax = plt.subplots(figsize=(6, 3))

# Row positions (y-values)
y_all  = 2
y_neg  = 1
y_pos  = 0

ax.scatter(pc1_scores,            np.full(105, y_all),
           c=colors,              s=20, marker="o", linewidths=0)
ax.scatter(pc1_scores[er_neg_mask], np.full(er_neg_mask.sum(), y_neg),
           c="black",             s=20, marker="o", linewidths=0)
ax.scatter(pc1_scores[er_pos_mask], np.full(er_pos_mask.sum(), y_pos),
           c="red",               s=20, marker="o", linewidths=0)

ax.set_yticks([0, 1, 2])
ax.set_yticklabels(["ER⁺", "ER⁻", "All"], fontsize=12)
ax.set_xlabel("Projection onto PC1", fontsize=12)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)

plt.tight_layout()
plt.savefig("figure1c_pc1_projection.png", dpi=150)
plt.close()
print("Saved → figure1c_pc1_projection.png")

print("\nAll figures saved successfully.")
