import zarr

store_path = "/home/vivanco/mif_resources/tests/data/microscopy/A1_F1.zarr"
root = zarr.open(store_path, mode='r')

# Print the visual hierarchy
print("Keys at root:", list(root.keys()))
print("\nFull Tree:")
print(root.tree())
