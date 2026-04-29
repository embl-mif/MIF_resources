from ome_zarr.io import parse_url
from ome_zarr.reader import Reader

path = "/home/vivanco/mif_resources/tests/data/microscopy/A1_F1.zarr"
reader = Reader(parse_url(path))
nodes = list(reader())

print(f"Found {len(nodes)} nodes.")
for n in nodes:
    print(f"Node path: {n.path} | Data shape: {n.data[0].shape}")


