# 📦 Import packages
import scanpy as sc
import anndata as ad
import celltypist
from celltypist import models
import mygene
import numpy as np
from scipy import sparse

# 📂 Load your dataset
adata = sc.read_h5ad("Chrishna-dataset2.h5ad")

# 🔍 Convert Ensembl IDs to gene symbols using mygene
mg = mygene.MyGeneInfo()
query = mg.querymany(adata.var_names.tolist(), scopes="ensembl.gene", fields="symbol", species="human")
mapping = {item['query']: item.get('symbol', item['query']) for item in query}
adata.var_names = [mapping.get(gene, gene) for gene in adata.var_names]

# 🧼 Normalize and log-transform the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 🧪 Save to .raw for CellTypist
adata.raw = adata.copy()

# 🧹 Clean NaNs in .raw.X
if sparse.issparse(adata.raw.X):
    raw_matrix = adata.raw.X.toarray()
else:
    raw_matrix = adata.raw.X
raw_matrix = np.nan_to_num(raw_matrix, nan=0.0)
adata.raw._X = raw_matrix

# 🧹 Clean NaNs in adata.X for PCA and UMAP
if sparse.issparse(adata.X):
    adata.X = adata.X.toarray()
adata.X = np.nan_to_num(adata.X, nan=0.0)

# 📉 Dimensionality reduction
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# 📥 Download and load CellTypist model
models.download_models()
model = models.Model.load(model="Immune_All_Low.pkl")

# 🧬 Run CellTypist annotation
results = celltypist.annotate(adata, model=model)

