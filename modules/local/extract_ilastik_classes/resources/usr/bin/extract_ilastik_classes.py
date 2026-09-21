#!/usr/bin/env python3
"""Read an ilastik .ilp project and print its class count/names as JSON.
Run upstream of PREALLOCATE_ZARR to resolve --num-channels/--channel-labels.
"""

import argparse
import json

from ilp_utils import extract_ilp_classes


def main():
    parser = argparse.ArgumentParser(description="Extract class/label info from an ilastik .ilp project.")
    parser.add_argument("--project", required=True, help="Path to Ilastik project file (.ilp)")
    args = parser.parse_args()

    classes = extract_ilp_classes(args.project)
    print(json.dumps({"num_channels": len(classes), "channel_labels": classes}))


if __name__ == "__main__":
    main()
