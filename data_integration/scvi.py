import numpy as np
import pandas as pd
import scanpy as sc
import scvi

# Load data
adata_merge = sc.read_h5ad("raw_rna.h5ad")
adata_merge.layers["counts"] = adata_merge.X.copy()

# Run scvi model
scvi.model.SCVI.setup_anndata(adata_merge, layer="counts", batch_key="sample")
model = scvi.model.SCVI(adata_merge, n_latent=30)
model.train(max_epochs = 100)

# Save results
scvi_embedding = model.get_latent_representation()
scvi_embedding = pd.DataFrame(scvi_embedding)
scvi_embedding.index = adata_merge.obs.index
scvi_embedding.to_csv("scvi.csv",index=True)
