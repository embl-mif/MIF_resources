#!/usr/bin/env python3

"""Concatenate every tile's `local_graph.json` contribution for one source
image and solve one multicut problem, producing a single global label
mapping (Module 5: LABEL_SOLVER).

Node IDs in every `local_graph.json` are already globally unique (offset by
`tile_index * label_offset_step` in Module 4), so no rejoining against tile
metadata is needed here -- this module just merges edges/disaffinities and
solves.
"""

import argparse
import json

import numpy as np

from label_matching import solve_global_multicut


def _load_local_graph(path: str) -> tuple[np.ndarray, np.ndarray]:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    # A tile with no edges round-trips its empty `edges` list as shape (0,),
    # not (0, 2) -- reshape explicitly so concatenation across tiles works.
    edges = np.asarray(data["edges"], dtype=np.uint64).reshape(-1, 2)
    disaffinities = np.asarray(data["disaffinities"], dtype=np.float32)
    return edges, disaffinities


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument("--local-graphs", type=str, nargs="+", required=True)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--output", type=str, required=True)

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    edge_list = []
    disaffinity_list = []
    for path in args.local_graphs:
        edges, disaffinities = _load_local_graph(path)
        edge_list.append(edges)
        disaffinity_list.append(disaffinities)

    label_mapping = solve_global_multicut(edge_list, disaffinity_list, beta=args.beta)

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(label_mapping, handle)

    print(f"Wrote {len(label_mapping)} label mapping(s) to {args.output}")


if __name__ == "__main__":
    main()

