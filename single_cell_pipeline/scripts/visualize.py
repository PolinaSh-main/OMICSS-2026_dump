import scanpy as sc

adata = sc.read_h5ad(
    "../results/harmony/harmony.h5ad"
)

print(adata)
print()
print("obs:")
print(adata.obs.columns)

print()
print("obsm:")
print(adata.obsm.keys())

print()
print("uns:")
print(adata.uns.keys())