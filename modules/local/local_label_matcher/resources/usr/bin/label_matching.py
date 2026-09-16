#!/usr/bin/env python3

"""Build one tile's contribution to a region-adjacency-graph (RAG) +
multicut stitching problem (Module 4: LOCAL_LABEL_MATCHER).

Graph node IDs are Cellpose label IDs offset by `tile_index * label_offset_step`
so they're globally unique across all tiles of one source image -- each
tile's Cellpose output independently restarts labeling at 1, and a boundary
edge connects a label from *this* tile to a label from a *neighbor* tile, so
both sides need to already agree on one global numbering before any edge is
built. That's why a neighbor's own `tile_index` has to be known here too, not
just this tile's -- offsetting a neighbor's labels with the wrong constant
would silently disagree with what the neighbor's own task uses when it later
contributes its own edges.

LABEL_SOLVER (Module 5) concatenates every tile's edges/disaffinities and
solves one multicut problem per source image.

The RAG-construction and multicut-solving core (`_halo_faces`,
`_add_boundary_edges`, `build_local_graph`, `_compute_edge_costs`,
`solve_global_multicut`, `relabel_segmentation`) is adapted from squirrel's
`stitching.py`
(https://github.com/embl-ic/squirrel/blob/main/squirrel/library/stitching.py,
GPL-3.0), with `_direction_axis` added to fix an axis-mapping bug in the
original that assumed 3D tiles unconditionally. squirrel is GPL-3.0 while
this repository is MIT-licensed -- flag this divergence to whoever owns
licensing decisions before distributing this file.
"""

import argparse
import json
from typing import Any

import numpy as np

LABEL_ARRAY_PATH = "scale0/nuclei_labels"

_DIRECTION_AXIS_3D = {"behind": 0, "bottom": 1, "right": 2}
_DIRECTION_AXIS_2D = {"bottom": 0, "right": 1}


def _direction_axis(direction: str, ndim: int) -> int:
    """Map a direction name to its array axis, for either a 2D or 3D tile.

    The original prototype hardcoded `{"behind": 0, "bottom": 1, "right": 2}`,
    which assumes `dims = (z, y, x)` -- correct for 3D tiles, but wrong for 2D
    ones (`dims = (y, x)`): `"right"` would map to axis 2, out of range, and
    `"bottom"` would map to axis 1, which for a 2D array is actually `x`.
    """
    if ndim == 3:
        mapping = _DIRECTION_AXIS_3D
    elif ndim == 2:
        mapping = _DIRECTION_AXIS_2D
    else:
        raise ValueError(f"Unsupported segmentation ndim={ndim}, expected 2 or 3")

    if direction not in mapping:
        raise ValueError(f"direction={direction!r} is not valid for a {ndim}D segmentation (no z axis)")

    return mapping[direction]


def _halo_faces(direction: str, ndim: int, overlap: tuple[int, ...]) -> tuple[tuple, tuple]:
    """
    Return matching halo slices for one neighbouring tile.

    Parameters
    ----------
    direction : {"right", "bottom", "behind"}
    ndim : int
    overlap : tuple[int]

    Returns
    -------
    this_face, neighbour_face : tuple[slice]
    """
    axis = _direction_axis(direction, ndim)

    this_face = []
    neighbour_face = []

    for ax in range(ndim):
        o = overlap[ax]

        if ax == axis:
            if o == 0:
                # No true overlap band to compare (adjacent tiles share zero
                # pixels), so fall back to a 1px touching-adjacency
                # comparison: this tile's last row/column against the
                # neighbour's first one. Weaker evidence than a real halo --
                # see label_matching.py's module docstring discussion -- but
                # strictly better than the empty-vs-empty no-op this used to
                # produce, which meant zero cross-tile edges were ever built
                # when overlap was 0. Also sidesteps that `slice(-0, None)`
                # is `slice(0, None)` (the whole axis), not an empty band,
                # since `-0 == 0` in Python.
                this_face.append(slice(-1, None))
                neighbour_face.append(slice(0, 1))
            else:
                # positive halo of this tile
                this_face.append(slice(-o, None))
                # negative halo of neighbour
                neighbour_face.append(slice(0, o))
        else:
            this_face.append(slice(None))
            neighbour_face.append(slice(None))

    return tuple(this_face), tuple(neighbour_face)


def _add_boundary_edges(edge_dict, seg, neighbour, direction, overlap, background):
    """Compute overlap-derived disaffinities for one neighbouring tile."""
    import bioimage_cpp as bic

    this_face, neighbour_face = _halo_faces(direction, seg.ndim, overlap)

    table = bic.utils.segmentation_overlap(
        np.ascontiguousarray(seg[this_face], dtype="uint64"),
        np.ascontiguousarray(neighbour[neighbour_face], dtype="uint64"),
    ).overlap_table()

    if table.shape[0] == 0:
        return

    # Compute normalized overlaps exactly as in bioimage-py

    label_a = table["label_a"]
    label_b = table["label_b"]
    counts = table["count"].astype(np.float64)

    # Total overlap size of each label_a
    unique_a, inverse = np.unique(label_a, return_inverse=True)

    size_a = np.zeros(len(unique_a), dtype=np.float64)
    np.add.at(size_a, inverse, counts)

    overlap_fraction = counts / size_a[inverse]

    # Convert to disaffinities

    for la, lb, frac in zip(label_a, label_b, overlap_fraction):
        if la == background or lb == background:
            continue

        edge = (min(int(la), int(lb)), max(int(la), int(lb)))

        # If multiple overlap entries exist for the same edge,
        # keep the strongest overlap (smallest disaffinity),
        # matching the original implementation.
        disaffinity = 1.0 - frac

        if edge in edge_dict:
            edge_dict[edge] = min(edge_dict[edge], disaffinity)
        else:
            edge_dict[edge] = disaffinity


def build_local_graph(
    seg,
    overlap,
    right=None,
    bottom=None,
    behind=None,
    background=0,
    default_disaffinity=0.9,
    n_workers=1,
):
    """
    Build the graph contribution of one tile.

    Parameters
    ----------
    seg : ndarray
        Tile segmentation including overlap, with node IDs already offset to
        be globally unique across the whole source image (see `_offset_labels`).
    overlap : tuple[int]
        Overlap width in pixels, in the same axis order as `seg.shape`
        (`(z, y, x)` for a 3D tile, `(y, x)` for a 2D one).
    right, bottom, behind : ndarray | None
        Forward neighbouring tiles, already offset the same way.

    Returns
    -------
    edges : (N, 2) ndarray
    disaffinities : (N,) ndarray
    """
    import bioimage_cpp as bic

    # Graph edges inside one tile get the default disaffinity

    rag = bic.graph.region_adjacency_graph(seg, number_of_threads=n_workers)

    edge_dict = {}

    for u, v in rag.uv_ids():
        if u == background or v == background:
            continue
        edge_dict[(min(u, v), max(u, v))] = default_disaffinity

    # Boundary edges get disaffinities based on overlap with objects in the neighbouring tile

    for direction, neighbour in (("right", right), ("bottom", bottom), ("behind", behind)):
        if neighbour is None:
            continue

        _add_boundary_edges(
            edge_dict=edge_dict,
            seg=seg,
            neighbour=neighbour,
            direction=direction,
            overlap=overlap,
            background=background,
        )

    edges = np.asarray(list(edge_dict.keys()), dtype=np.uint64)
    if len(edges) == 0:
        edges = np.empty((0, 2), dtype=np.uint64)

    disaffinities = np.asarray(list(edge_dict.values()), dtype=np.float32)

    return edges, disaffinities


def _compute_edge_costs(disaffinities: np.ndarray, beta: float = 0.5, eps: float = 1e-6) -> np.ndarray:
    """Convert disaffinities in [0, 1] into multicut edge costs.

    Positive cost => attractive (merge), negative => repulsive (keep separate).
    `beta` adds a uniform shift independent of the observed disaffinity: at
    beta=0.5 it contributes nothing; raising beta above 0.5 shifts every edge's
    cost down (towards repulsive), biasing the solve towards more, smaller
    separate objects; lowering it below 0.5 shifts costs up, biasing towards
    merging more tiles' labels together. Verified empirically (not just
    derived), since this formula is unchanged from the original prototype.
    """
    p = np.clip(disaffinities, eps, 1.0 - eps)

    costs = np.log((1.0 - p) / p)
    costs += np.log((1.0 - beta) / beta)

    return costs.astype(np.float32)


def solve_global_multicut(
    edge_list: list[np.ndarray],
    disaffinity_list: list[np.ndarray],
    beta: float = 0.5,
) -> dict[int, int]:
    """Solve one multicut problem over every tile's local graph contribution
    (Module 5: LABEL_SOLVER).

    Parameters
    ----------
    edge_list : list[np.ndarray]
        Per-tile edge arrays from `build_local_graph`, each shape (N, 2), node
        IDs already globally unique (offset by `tile_index * label_offset_step`).
    disaffinity_list : list[np.ndarray]
        Corresponding per-tile disaffinity arrays.
    beta : float
        Prior parameter for the multicut cost transform.

    Returns
    -------
    label_mapping : dict[int, int]
        Global node ID -> global label, for every node ID that appears in at
        least one edge. A node that never touches another label (no in-tile
        neighbor, no boundary overlap) has no entry here -- callers must fall
        back to the node's own ID unchanged (it's already globally unique via
        the offset scheme and needs no merging), not treat a missing key as
        an error. See Module 6's `relabel_segmentation`.
    """
    import nifty.graph as ng
    import nifty.graph.opt.multicut as nmc

    # Merge all local graphs
    edges = np.concatenate(edge_list, axis=0)

    if len(edges) == 0:
        return {}

    disaffinities = np.concatenate(disaffinity_list)

    # Remove duplicate edges (a boundary edge is built independently by both
    # tiles on either side of it, so it appears twice). Keep the strongest
    # attraction (= smallest disaffinity).
    edge_dict: dict[tuple[int, int], float] = {}
    for edge, dis in zip(edges, disaffinities):
        key = (min(int(edge[0]), int(edge[1])), max(int(edge[0]), int(edge[1])))

        if key in edge_dict:
            edge_dict[key] = min(edge_dict[key], dis)
        else:
            edge_dict[key] = dis

    # Convert back to arrays
    edges = np.asarray(list(edge_dict.keys()), dtype=np.uint64)
    disaffinities = np.asarray(list(edge_dict.values()), dtype=np.float32)

    # Relabel graph nodes to a dense index set
    graph_labels = np.unique(edges.ravel())
    dense_edges = np.searchsorted(graph_labels, edges)

    # Build graph
    graph = ng.UndirectedGraph(len(graph_labels))
    graph.insertEdges(dense_edges)

    # Compute multicut costs
    costs = _compute_edge_costs(disaffinities, beta=beta)

    # Solve
    objective = nmc.multicutObjective(graph, costs)
    solver = objective.kernighanLinFactory().create(objective)
    node_labels = solver.optimize()

    # Compress component ids to consecutive labels
    unique_components = np.unique(node_labels)
    component_to_label = {component: i + 1 for i, component in enumerate(unique_components)}

    # Final mapping: global node ID -> stitched global label
    label_mapping = {
        int(label): component_to_label[component]
        for label, component in zip(graph_labels, node_labels)
    }

    return label_mapping


def relabel_segmentation(
    seg: np.ndarray,
    label_mapping: dict[int, int],
    background: int = 0,
) -> np.ndarray:
    """Apply a global label mapping to one tile's segmentation (Module 6b:
    WRITE_RELABELED_TILE).

    `seg` must already be in global-node-id space (i.e. already passed
    through `_offset_labels` with this tile's own `tile_index`) -- the
    `label_mapping.get(int(label), int(label))` fallback below is only
    correct because a node ID absent from `label_mapping` is, by
    construction, already globally unique (see Module 5's
    `solve_global_multicut` docstring), not because any arbitrary local
    label happens to be globally unique on its own.

    Parameters
    ----------
    seg : np.ndarray
    label_mapping : dict[int, int]
        Output of `solve_global_multicut`.
    background : int

    Returns
    -------
    relabeled : np.ndarray
    """
    labels, inverse = np.unique(seg, return_inverse=True)

    new_labels = labels.copy()

    for i, label in enumerate(labels):
        if label == background:
            continue
        new_labels[i] = label_mapping.get(int(label), int(label))

    return new_labels[inverse].reshape(seg.shape)


def _offset_labels(seg: np.ndarray, offset: int, background: int = 0) -> np.ndarray:
    """Shift every non-background label by `offset` to make it globally unique.

    Kept as uint32 to match what the original prototype's manual `+ 2**16`
    offset test exercised against `bic.graph.region_adjacency_graph` -- callers
    must keep `label_offset_step * max(tile_index) + max(local_label)` under
    2**32, which is why this raises rather than silently wrapping.
    """
    max_value = (int(seg.max()) if seg.size else 0) + offset
    if max_value > np.iinfo(np.uint32).max:
        raise ValueError(
            f"Offset label {max_value} exceeds uint32 range ({np.iinfo(np.uint32).max}) -- "
            "reduce --label-offset-step or the number of tiles per source image"
        )

    offset_seg = seg.astype(np.uint32)
    mask = offset_seg != background
    offset_seg[mask] = offset_seg[mask] + np.uint32(offset)
    return offset_seg


def _open_label_array(path: str):
    import zarr
    array = zarr.open(path, mode="r")[LABEL_ARRAY_PATH]
    if array.ndim > 3 and array.shape[0] == 1:
        # Singleton leading axis from an OME-Zarr `t` (or stray) dimension the
        # tiled-stitching pipeline doesn't model -- drop it here so every
        # downstream ndim==2/3 check sees the real spatial rank. A genuine
        # multi-frame axis (shape[0] > 1) is left untouched and will still
        # correctly raise, since this pipeline has no per-timepoint stitching.
        array = array[0]
    return array


def _load_full(path: str) -> np.ndarray:
    return np.asarray(_open_label_array(path)[:])


def _load_halo_band(path: str, direction: str, ndim: int, overlap: tuple[int, ...]) -> np.ndarray:
    """Read only the shared halo band from a neighbor's zarr store -- never
    the full neighbor tile. Safe to re-slice inside `_add_boundary_edges`
    (it recomputes the same `neighbour_face` and reapplies it): this band's
    size along the tested axis already equals `overlap[axis]` (or 1, for the
    touching-adjacency fallback `_halo_faces` uses when `overlap[axis]` is
    0), so re-slicing with the same face on it is a no-op."""
    _, neighbour_face = _halo_faces(direction, ndim, overlap)
    array = _open_label_array(path)
    return np.ascontiguousarray(array[neighbour_face])


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument("--self-path", type=str, required=True)
    parser.add_argument("--tile-index", type=int, required=True)
    parser.add_argument("--label-offset-step", type=int, required=True)

    parser.add_argument("--tile-overlap-x", type=int, required=True)
    parser.add_argument("--tile-overlap-y", type=int, required=True)
    parser.add_argument("--tile-overlap-z", type=int, default=None)

    parser.add_argument("--right-path", type=str, default=None)
    parser.add_argument("--right-tile-index", type=int, default=None)
    parser.add_argument("--bottom-path", type=str, default=None)
    parser.add_argument("--bottom-tile-index", type=int, default=None)
    parser.add_argument("--behind-path", type=str, default=None)
    parser.add_argument("--behind-tile-index", type=int, default=None)

    parser.add_argument("--background", type=int, default=0)
    parser.add_argument("--default-disaffinity", type=float, default=0.9)

    parser.add_argument("--output", type=str, required=True)

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    self_seg = _load_full(args.self_path)
    ndim = self_seg.ndim

    if ndim == 3:
        if args.tile_overlap_z is None:
            raise ValueError("--tile-overlap-z is required for a 3D segmentation")
        overlap = (args.tile_overlap_z, args.tile_overlap_y, args.tile_overlap_x)
    elif ndim == 2:
        overlap = (args.tile_overlap_y, args.tile_overlap_x)
    else:
        raise ValueError(f"Unsupported segmentation ndim={ndim}, expected 2 or 3")

    self_seg = _offset_labels(self_seg, args.tile_index * args.label_offset_step, args.background)

    neighbor_args = {
        "right": (args.right_path, args.right_tile_index),
        "bottom": (args.bottom_path, args.bottom_tile_index),
        "behind": (args.behind_path, args.behind_tile_index),
    }

    neighbors: dict[str, Any] = {}
    for direction, (path, neighbor_tile_index) in neighbor_args.items():
        if path is None:
            continue
        if neighbor_tile_index is None:
            raise ValueError(f"--{direction}-path was given without --{direction}-tile-index")

        band = _load_halo_band(path, direction, ndim, overlap)
        neighbors[direction] = _offset_labels(
            band, neighbor_tile_index * args.label_offset_step, args.background
        )

    edges, disaffinities = build_local_graph(
        seg=self_seg,
        overlap=overlap,
        right=neighbors.get("right"),
        bottom=neighbors.get("bottom"),
        behind=neighbors.get("behind"),
        background=args.background,
        default_disaffinity=args.default_disaffinity,
    )

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump({"edges": edges.tolist(), "disaffinities": disaffinities.tolist()}, handle)

    print(f"Wrote {len(edges)} edge(s) to {args.output}")


if __name__ == "__main__":
    main()

