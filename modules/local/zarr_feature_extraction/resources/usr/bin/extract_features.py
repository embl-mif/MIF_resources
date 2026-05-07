#!/usr/bin/env python3
import sys
import os

lib_path = os.environ.get('MODULE_LIB')

# Add the lib path so Python can see the 'zarr_features' folder
sys.path.append(lib_path)

from zarr_features.cli import main

if __name__ == "__main__":
    sys.exit(main())
