import zarr
import os

base_path = "/home/vivanco/mif_resources/tests/data/microscopy/A1_F1.zarr"
label_dir = os.path.join(base_path, "labels")

print(f"DEBUG: Checking path: {label_dir}")
print(f"DEBUG: Does directory exist? {os.path.exists(label_dir)}")
print(f"DEBUG: Does zarr.json exist inside? {os.path.exists(os.path.join(label_dir, 'zarr.json'))}")

if os.path.exists(os.path.join(label_dir, 'zarr.json')):
    # Try reading the metadata to see if it's a valid zarr group
    with open(os.path.join(label_dir, 'zarr.json'), 'r') as f:
        print("DEBUG: zarr.json content preview:", f.read(100))
