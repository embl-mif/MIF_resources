from ome_zarr.io import parse_url
from ome_zarr.reader import Reader

path = '/home/vivanco/mif_resources/tests/data/microscopy/A1_F1.zarr'
reader = Reader(parse_url(path))

for i, node in enumerate(reader()):
    print(f"--- Node {i} ---")
    # This prints all attributes and methods available on the node
    print(dir(node))
