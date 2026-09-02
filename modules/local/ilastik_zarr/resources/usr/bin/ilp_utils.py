"""Shared helpers for reading Ilastik .ilp project files."""

import math
import os

import numpy as np
import h5py


def extract_ilp_classes(ilp_path: str):
    """get class names out of the .ilp file. number of classes/labels = number of output channels."""
    if not os.path.exists(ilp_path):
        raise FileNotFoundError(f"Ilastik project file not found: {ilp_path}")

    with h5py.File(ilp_path, "r") as f:
        if "PixelClassification/LabelNames" not in f:
            raise KeyError("Could not find '/PixelClassification/LabelNames' in the .ilp file.")

        raw_labels = f["PixelClassification/LabelNames"][()]
        return [l.decode("utf-8") if isinstance(l, bytes) else str(l) for l in raw_labels]


def _ilp_max_sigma(f: h5py.File) -> float:
    """biggest sigma across the enabled feature scales. falls back to 1.6 if nothing's set."""
    if "FeatureSelections/Scales" not in f or "FeatureSelections/SelectionMatrix" not in f:
        return 1.6

    scales = np.asarray(f["FeatureSelections/Scales"][()], dtype=float).ravel()
    matrix = np.asarray(f["FeatureSelections/SelectionMatrix"][()])

    if matrix.size == 0 or not np.any(matrix):
        return 1.6

    active_indices = np.where(matrix.any(axis=0))[0]
    if len(active_indices) == 0:
        return 1.6

    return float(np.max(scales[active_indices]))


def inspect_ilp(project_path: str):
    """check if the project is 2D or 3D, and work out a halo size from the feature sigmas."""
    with h5py.File(project_path, "r") as f:
        is_2d = False
        if "FeatureSelections/ComputeIn2d" in f:
            val = f["FeatureSelections/ComputeIn2d"][()]
            is_2d = bool(np.any(val))

        max_sigma = _ilp_max_sigma(f)
        recommended_halo = int(math.ceil(3.5 * max_sigma))

    return is_2d, recommended_halo
