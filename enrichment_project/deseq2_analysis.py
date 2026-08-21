import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt

# ============================================================
# 1. Загрузка данных
# ============================================================
counts = pd.read_csv("dataset/STAR_counts.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("dataset/metadata.tsv", sep="\t", index_col=0)
gene_ann = pd.read_csv("dataset/gencode.v42.genes.tsv", sep="\t", index_col=0)

# pydeseq2 хочет samples x genes (транспонированную матрицу)
counts_T = counts.T

# ============================================================
# 2. Pre-filtering (как в R: оставляем гены с >= 10 counts в >= 3 samples)
# ============================================================
keep = (counts_T >= 10).sum(axis=0) >= 3
counts_T = counts_T.loc[:, keep]

# ============================================================
# 3. DESeqDataSet + нормализация + DESeq
# ============================================================
dds = DeseqDataSet(
    counts=counts_T,
    metadata=metadata,
    design="~Subtype",
)

dds.fit_size_factors()
print("Size factors:")
print(dds.obs["size_factors"])

# Нормализованные counts (для боксплотов и PCA)
norm_counts = counts_T / dds.obs["size_factors"].values[:, None]

# ============================================================
# 4. Боксплот ESR1
# ============================================================
gene_name = "ESR1"
gene_id = gene_ann[gene_ann["gene_name"] == gene_name].index[0]

if gene_id in norm_counts.columns:
    plot_df = pd.DataFrame({
        "counts": norm_counts[gene_id],
        "Subtype": metadata.loc[norm_counts.index, "Subtype"],
    })
    plot_df.boxplot(column="counts", by="Subtype", figsize=(8, 5))
    plt.title("ESR1")
    plt.suptitle("")
    plt.ylabel("Normalized counts")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("esr1_boxplot.png", dpi=150)
    plt.show()

# ============================================================
# 5. Differential expression: Luminal A vs TNBC
# ============================================================
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["Subtype", "LuminalA", "TNBC"])
stat_res.summary()

res_raw = stat_res.results_df.copy()

# Shrinkage - пока заглушка
res_shrink = res_raw.copy()

# ============================================================
# 6. Добавляем аннотацию генов
# ============================================================
res_raw = res_raw.join(gene_ann, how="left")
res_shrink = res_shrink.join(gene_ann, how="left")

# Сортируем по padj
res_raw.sort_values("padj", inplace=True)
res_shrink.sort_values("padj", inplace=True)

# ============================================================
# 7. Боксплоты для top-4 DEGs
# ============================================================
top4 = res_raw.head(4)
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for ax, (gene_id, row) in zip(axes, top4.iterrows()):
    if gene_id in norm_counts.columns:
        plot_df = pd.DataFrame({
            "counts": norm_counts[gene_id],
            "Subtype": metadata.loc[norm_counts.index, "Subtype"],
        })
        plot_df.boxplot(column="counts", by="Subtype", ax=ax)
        ax.set_title(row.get("gene_name", gene_id))
        ax.set_xlabel("")
        plt.sca(ax)
        plt.xticks(rotation=45, fontsize=7)
    ax.set_ylabel("")
plt.suptitle("Top 4 DEGs by padj")
plt.tight_layout()
plt.savefig("top4_degs.png", dpi=150)
plt.show()

# ============================================================
# 8. Volcano plot (простой, без EnhancedVolcano)
# ============================================================
def volcano(df, title="Volcano", fc_cut=1, p_cut=0.05):
    df = df.dropna(subset=["padj", "log2FoldChange"])
    neg_log_p = -np.log10(df["padj"].clip(lower=1e-300))
    colors = np.where(
        (df["padj"] < p_cut) & (df["log2FoldChange"].abs() >= fc_cut),
        np.where(df["log2FoldChange"] > 0, "firebrick", "dodgerblue"),
        "grey"
    )
    plt.figure(figsize=(8, 6))
    plt.scatter(df["log2FoldChange"], neg_log_p, c=colors, s=5, alpha=0.6)
    plt.axhline(-np.log10(p_cut), color="grey", linestyle="--", linewidth=0.5)
    plt.axvline(-fc_cut, color="grey", linestyle="--", linewidth=0.5)
    plt.axvline(fc_cut, color="grey", linestyle="--", linewidth=0.5)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(padj)")
    plt.title(title)
    plt.xlim(-10, 10)
    plt.tight_layout()
    plt.savefig(f"{title.replace(' ', '_').lower()}.png", dpi=150)
    plt.show()

volcano(res_raw, "Volcano (raw)")
volcano(res_shrink, "Volcano (shrunk)")

# ============================================================
# 9. Фильтрация значимых генов
# ============================================================
res_signif = res_shrink.dropna(subset=["padj"])
res_signif = res_signif[
    (res_signif["padj"] < 0.05) &
    (res_signif["log2FoldChange"].abs() >= 1) &
    (res_signif["baseMean"] >= 100) &
    (res_signif["gene_type"] == "protein_coding")
]
res_signif_up = res_signif[res_signif["log2FoldChange"] > 0]
res_signif_down = res_signif[res_signif["log2FoldChange"] < 0]

res_signif_up.to_csv("overexpressed_genes.csv")
res_signif_down.to_csv("downexpressed_genes.csv")

print(f"Up: {len(res_signif_up)}, Down: {len(res_signif_down)}")

# ============================================================
# 10. Assignment 1: scatter raw FC vs shrunk FC, colored by baseMean
# ============================================================
merged = res_raw[["log2FoldChange", "baseMean"]].rename(
    columns={"log2FoldChange": "lfc_raw"}
).join(
    res_shrink[["log2FoldChange"]].rename(columns={"log2FoldChange": "lfc_shrunk"})
).dropna()

plt.figure(figsize=(7, 6))
plt.scatter(merged["lfc_raw"], merged["lfc_shrunk"],
            c=np.log10(merged["baseMean"] + 1), s=3, alpha=0.5, cmap="viridis")
plt.colorbar(label="log10(baseMean)")
plt.xlabel("Raw log2FC")
plt.ylabel("Shrunken log2FC")
plt.title("Raw vs Shrunken LFC")
plt.tight_layout()
plt.savefig("raw_vs_shrunk_lfc.png", dpi=150)
plt.show()

res_raw.to_csv("deseq2_full_results.csv")