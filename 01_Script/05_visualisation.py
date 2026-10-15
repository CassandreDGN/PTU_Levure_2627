import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from pathlib import Path
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

# --- 1. CONFIGURATION DES CHEMINS ---
MATRIX_FILE = Path("/data/projet1/02_Pipeline_Output/6_Diploid_Genotypes/diploid_matrix_all_strains.tsv")
OUT_DIR = Path("/data/projet1/02_Pipeline_Output/7_Synthesis_and_Plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=== ÉTAPE 5 : SYNTHÈSE DES DONNÉES ET VISUALISATIONS (COLORBLIND FRIENDLY) ===")

if not MATRIX_FILE.exists():
    raise FileNotFoundError(f"Fichier de matrice introuvable : {MATRIX_FILE}")

print("Chargement de la matrice diploïde globale...")
df = pd.read_csv(MATRIX_FILE, sep="\t", index_col="protein")
n_genes, n_strains = df.shape
print(f"-> Matrice chargée : {n_genes} gènes x {n_strains} souches.")

# --- 2. BILANS CHIFFRÉS ---
print("\\n1/4. Calcul des bilans statistiques...")
strain_stats = []
for strain in df.columns:
    col = df[strain]
    strain_stats.append({
        "strain": strain,
        "n_homo_functional_2": (col == 2).sum(),
        "n_hetero_1": (col == 1).sum(),
        "n_inactive_loss_0": (col == 0).sum(),
        "pct_inactive": round(((col == 0).sum() / n_genes) * 100, 2)
    })
df_strain_summary = pd.DataFrame(strain_stats)
df_strain_summary.to_csv(OUT_DIR / "strain_summary.tsv", sep="\t", index=False)

gene_stats = []
for gene, row in df.iterrows():
    gene_stats.append({
        "protein": gene,
        "strains_homo_2": (row == 2).sum(),
        "strains_hetero_1": (row == 1).sum(),
        "strains_inactive_0": (row == 0).sum(),
        "pct_strains_inactive": round(((row == 0).sum() / n_strains) * 100, 2)
    })
df_gene_summary = pd.DataFrame(gene_stats)
df_gene_summary.to_csv(OUT_DIR / "gene_summary_global.tsv", sep="\t", index=False)

# --- 3. CONFIGURATION VISUELLE (STYLE ACADÉMIQUE COLORBLIND) ---
sns.set_theme(style="ticks", context="paper", font_scale=1.2)

# Palette Okabe-Ito : standard scientifique pour l'accessibilité visuelle
COLOR_LOSS = "#D55E00"   # Vermillon
COLOR_HETERO = "#E69F00" # Orange
COLOR_FUNC = "#0072B2"   # Bleu

print("\\n2/4. Génération des graphiques (PNG + SVG)...")

# A. Heatmap Hiérarchique (Clustermap)
df_variants = df[(df != 2).any(axis=1)]
n_var_genes = len(df_variants)
print(f"  -> Heatmap ({n_var_genes} gènes variants)...")

if n_var_genes > 0:
    cat_cmap = ListedColormap([COLOR_LOSS, COLOR_HETERO, COLOR_FUNC])
    
    g = sns.clustermap(
        df_variants,
        cmap=cat_cmap,
        linewidths=0,
        figsize=(10, 8),
        cbar_pos=None,
        yticklabels=False,
        xticklabels=True,
        tree_kws={"linewidths": 1.2}
    )
    g.fig.suptitle(f"Profils d'inactivation diploïde ({n_var_genes} gènes variants)", y=1.05, fontweight='bold')
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=9)
    
    legend_patches = [
        mpatches.Patch(color=COLOR_LOSS, label='Perte (0)'),
        mpatches.Patch(color=COLOR_HETERO, label='Hétérozygote (1)'),
        mpatches.Patch(color=COLOR_FUNC, label='Fonctionnel (2)')
    ]
    g.fig.legend(handles=legend_patches, loc='center right', bbox_to_anchor=(1.15, 0.5), title="Génotype", frameon=False)
    
    g.savefig(OUT_DIR / "01_heatmap_clustermap.png", dpi=300, bbox_inches="tight")
    g.savefig(OUT_DIR / "01_heatmap_clustermap.svg", format="svg", bbox_inches="tight")
    plt.close()

# B. PCA Enrichie avec Charge Mutationnelle (Palette séquentielle 'viridis')
print("  -> PCA spatiale enrichie...")
pca = PCA(n_components=2)
pca_res = pca.fit_transform(df.T.values)

df_pca = pd.DataFrame({"PC1": pca_res[:, 0], "PC2": pca_res[:, 1], "strain": df.columns})
df_pca = pd.merge(df_pca, df_strain_summary, on="strain")

var_exp = pca.explained_variance_ratio_ * 100

fig, ax = plt.subplots(figsize=(9, 7))
sc = sns.scatterplot(
    data=df_pca, x="PC1", y="PC2", 
    hue="n_inactive_loss_0", palette="viridis", 
    size="n_inactive_loss_0", sizes=(80, 250), 
    edgecolor="black", zorder=3, ax=ax
)

for _, row in df_pca.iterrows():
    ax.annotate(row["strain"], (row["PC1"], row["PC2"]), xytext=(7, 7), textcoords="offset points", fontsize=9)

ax.set_title("Espace Génomique des Souches (PCA)", fontweight='bold', pad=15)
ax.set_xlabel(f"Composante 1 ({var_exp[0]:.1f}% variance)")
ax.set_ylabel(f"Composante 2 ({var_exp[1]:.1f}% variance)")

sc.legend_.set_title("Gènes Inactifs")
sns.despine()
plt.tight_layout()

plt.savefig(OUT_DIR / "02_pca_strains.png", dpi=300, bbox_inches="tight")
plt.savefig(OUT_DIR / "02_pca_strains.svg", format="svg", bbox_inches="tight")
plt.close()

# C. Barplot Empilé (Charge Mutationnelle Totale)
print("  -> Stacked Barplot des altérations...")
df_sorted = df_strain_summary.sort_values(by=["n_inactive_loss_0", "n_hetero_1"], ascending=[False, False])
df_sorted.set_index("strain", inplace=True)

fig, ax = plt.subplots(figsize=(12, 6))
df_sorted[['n_hetero_1', 'n_inactive_loss_0']].plot(
    kind='bar', stacked=True, 
    color=[COLOR_HETERO, COLOR_LOSS], 
    width=0.8, ax=ax, edgecolor="black", linewidth=0.5
)

ax.set_title("Charge mutationnelle par souche (Allèles altérés)", fontweight='bold', pad=15)
ax.set_xlabel("Souche")
ax.set_ylabel("Nombre de gènes")
ax.legend(["Hétérozygote (1)", "Perte de fonction (0)"], frameon=False)
plt.xticks(rotation=45, ha="right")
sns.despine()

for container in ax.containers:
    if container == ax.containers[-1]:
        ax.bar_label(container, padding=3, fmt='%.0f', fontsize=9, color=COLOR_LOSS, fontweight='bold')

plt.tight_layout()
plt.savefig(OUT_DIR / "03_strain_mutational_burden.png", dpi=300, bbox_inches="tight")
plt.savefig(OUT_DIR / "03_strain_mutational_burden.svg", format="svg", bbox_inches="tight")
plt.close()

# D. Top Gènes Impactés (Horizontal avec valeurs)
print("  -> Barplot des gènes d'intérêt...")
top_genes = df_gene_summary.sort_values(by="strains_inactive_0", ascending=False).head(20)

fig, ax = plt.subplots(figsize=(10, 7))
sns.barplot(data=top_genes, x="strains_inactive_0", y="protein", color=COLOR_FUNC, edgecolor="black", ax=ax)

ax.bar_label(ax.containers[0], padding=5, fontsize=10)

ax.set_title("Top 20 des gènes les plus inactivés dans la population", fontweight='bold', pad=15)
ax.set_xlabel("Nombre de souches présentant une perte totale (Score = 0)")
ax.set_ylabel("Protéine / Gène")
sns.despine()
plt.tight_layout()

plt.savefig(OUT_DIR / "04_top_inactive_genes.png", dpi=300, bbox_inches="tight")
plt.savefig(OUT_DIR / "04_top_inactive_genes.svg", format="svg", bbox_inches="tight")
plt.close()

print("\\n✅ Visualisations de qualité académique générées avec succès !")