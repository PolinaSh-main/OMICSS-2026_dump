import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# ============================================================
# 1. Загрузка — напрямую из результатов первого скрипта
# ============================================================
# Если запускаешь отдельным файлом, просто читай CSV,
# которые сохранил первый скрипт:
res_signif_up = pd.read_csv("overexpressed_genes.csv", index_col=0)
res_signif_down = pd.read_csv("downexpressed_genes.csv", index_col=0)

# Для GSEA нужен полный (нефильтрованный) результат —
# либо передай res_raw напрямую, либо сохрани его в первом скрипте:
#   res_raw.to_csv("deseq2_full_results.csv")
res_raw = pd.read_csv("deseq2_full_results.csv", index_col=0)

# Для PCA — нормализованные counts и metadata
# Добавь в первый скрипт:
#   norm_counts.to_csv("norm_counts.csv")
counts = pd.read_csv("dataset/STAR_counts.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("dataset/metadata.tsv", sep="\t", index_col=0)

# ============================================================
# 2. ORA
# ============================================================
dbs = [
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
]

ora_up = gp.enrichr(
    gene_list=res_signif_up["gene_name"].dropna().tolist(),
    gene_sets=dbs,
    organism="human",
    outdir=None,
)
ora_down = gp.enrichr(
    gene_list=res_signif_down["gene_name"].dropna().tolist(),
    gene_sets=dbs,
    organism="human",
    outdir=None,
)

ora_up_sig = ora_up.results[ora_up.results["Adjusted P-value"] < 0.05]
ora_down_sig = ora_down.results[ora_down.results["Adjusted P-value"] < 0.05]

print("=== ORA UP ===")
print(ora_up_sig.head(20).to_string())
print("\n=== ORA DOWN ===")
print(ora_down_sig.head(20).to_string())

# ============================================================
# 3. GSEA
# ============================================================
res_raw = res_raw.dropna(subset=["stat"])
res_raw.sort_values("stat", ascending=False, inplace=True)

rnk = res_raw[["gene_name", "stat"]].dropna()
rnk.columns = ["gene", "score"]
rnk = rnk.sort_values("score", ascending=False)

gsea_res = gp.prerank(
    rnk=rnk,
    gene_sets="gobp_symbols.gmt",
    min_size=15,
    max_size=500,
    permutation_num=100,
    seed=54321,
    outdir=None,
    threads=1,
)

gsea_df = gsea_res.res2d.copy()
gsea_df["NES"] = gsea_df["NES"].astype(float)
gsea_df["FDR q-val"] = gsea_df["FDR q-val"].astype(float)
gsea_signif = gsea_df[gsea_df["FDR q-val"] < 0.05].sort_values("NES", ascending=False)

total_up = (gsea_signif["NES"] > 0).sum()
total_down = (gsea_signif["NES"] < 0).sum()
print(f"\nSignificant: {len(gsea_signif)} (Up={total_up}, Down={total_down})")

# --- Summary plot ---
top = pd.concat([gsea_signif.head(10), gsea_signif.tail(10)])
top["Enrichment"] = np.where(top["NES"] > 0, "Up-regulated", "Down-regulated")

# размер точки — кол-во генов в пересечении
size_col = "matched_size" if "matched_size" in top.columns else "Tag %"
top["dot_size"] = pd.to_numeric(top[size_col], errors="coerce").fillna(30)

colors = {"Up-regulated": "firebrick", "Down-regulated": "dodgerblue"}
fig, ax = plt.subplots(figsize=(10, 8))
for label, grp in top.groupby("Enrichment"):
    ax.scatter(grp["NES"], grp["Term"], s=grp["dot_size"] * 3,
               color=colors[label], label=label, edgecolors="black", alpha=0.8)
ax.axvline(0, color="grey", linewidth=0.8)
ax.set_xlabel("Normalized Enrichment Score")
ax.set_title(f"Top 10 (Total: Up={total_up}, Down={total_down})")
ax.legend()
plt.tight_layout()
plt.savefig("gsea_summary.png", dpi=150)
plt.show()

# --- Enrichment plot для EMT ---
try:
    match = [t for t in gsea_signif["Term"] if "mesenchymal" in t.lower()]
    if match:
        from gseapy import gseaplot
        gseaplot(rank_metric=gsea_res.ranking,
                 term=match[0],
                 ofname="enrichment_emt.png",
                 **gsea_res.results[match[0]])
except Exception as e:
    print(f"EMT enrichment plot skipped: {e}")

# ============================================================
# 4. PCA
# ============================================================
# Берём сырые counts, делаем простую log-нормализацию
log_counts = np.log2(counts + 1)

pca = PCA(n_components=2)
pcs = pca.fit_transform(log_counts.T)

pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=counts.columns)
pca_df["Subtype"] = metadata.loc[pca_df.index, "Subtype"]

fig, ax = plt.subplots(figsize=(8, 6))
for subtype, grp in pca_df.groupby("Subtype"):
    ax.scatter(grp["PC1"], grp["PC2"], label=subtype, s=60, alpha=0.8)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA Plot")
ax.legend()
ax.grid(color="grey", alpha=0.3)
plt.tight_layout()
plt.savefig("pca_plot.png", dpi=150)
plt.show()