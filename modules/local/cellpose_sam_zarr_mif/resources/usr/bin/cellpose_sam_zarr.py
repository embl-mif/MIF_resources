#!/usr/bin/env python3

"""
Instructions to run this code:

* Setup environment for cellpose and ngff-zarr

Using pip:
source /g/cba/miniconda3/init_conda.sh 
conda create --prefix /g/cba/miniconda3/envs/cellpose_zarr -c conda-forge --override-channels python=3.12 -y
conda activate cellpose_zarr
pip install "ngff-zarr[all]" cellpose
pip install fsspec[http]


* Run code on EMBL cluster CLI:

srun --nodes=1 --cpus-per-task=4 --mem-per-cpu=32000 --time=01:00:00 -G 1 -C "gpu=3090" -p gpu-el8 --pty /bin/bash
conda activate cellpose_zarr
python bin/cellpose_sam_zarr.py --channels 4 -i /g/cba/exchange/ipf-tma-analysis-test/ome-zarr/Run01_Point0001_Crop.ome.zarr -o /g/cba/exchange/ipf-tma-analysis-test/ome-zarr/Run01_Point0001_Crop_Labels_Tischi_2.ome.zarr
python bin/cellpose_sam_zarr.py --channels 4 -i /g/cba/exchange/ipf-tma-analysis-test/ome-zarr/Run01_Point0001.ome.zarr -o /g/cba/exchange/ipf-tma-analysis-test/ome-zarr/Run01_Point0001_Labels_Tischi.ome.zarr

python bin/cellpose_sam_zarr.py -i data/xy-nuclei-v4.ome.zarr -o data/xy-nuclei-labels-v4.ome.zarr
python bin/cellpose_sam_zarr.py -i data/xyzct-v5.ome.zarr -o data/xyzct-labels-v5.ome.zarr
"""

import argparse
import math
import os
import re
from typing import Any

import ngff_zarr as nz
import numpy as np
import zarr
from cellpose import core, io, models
from ngff_zarr.methods import Methods


def normalize_ome_zarr_version(version: str) -> str:
    normalized = str(version).strip()
    if normalized in {"0.4", "0.5"}:
        return normalized
    raise ValueError(
        f"Unsupported OME-Zarr version '{version}'. Supported values are: 0.4, 0.5"
    )


def detect_ome_zarr_version(input_zarr: str) -> str:
    """Detect OME-Zarr version from root attributes for both v0.4 and v0.5 layouts."""
    z_group = zarr.open(input_zarr, mode="r")
    attrs = z_group.attrs

    ome = attrs.get("ome")
    if isinstance(ome, dict) and "version" in ome:
        return normalize_ome_zarr_version(ome["version"])

    multiscales = attrs.get("multiscales", [{}])
    if multiscales and isinstance(multiscales[0], dict) and "version" in multiscales[0]:
        return normalize_ome_zarr_version(multiscales[0]["version"])

    return "0.4"


def scale_factors_from_input(multiscales: Any, label_dims: list[str], base_level: int) -> list[Any]:
    images = getattr(multiscales, "images", None) or []
    if base_level < 0 or base_level >= len(images):
        return []

    if len(images) <= base_level + 1:
        return []

    base_image = images[base_level]
    base_scales = {dim: float(base_image.scale.get(dim, 1.0)) for dim in label_dims}
    inferred: list[dict[str, int]] = []
    for image in images[base_level + 1:]:
        factors_for_level: dict[str, int] = {}
        for dim in label_dims:
            # to_multiscales expects spatial scaling only.
            if dim in {"t", "c"}:
                continue
            base = base_scales.get(dim, 1.0)
            current = float(image.scale.get(dim, base))
            factor = int(round(current / base)) if base != 0.0 else 1
            factors_for_level[dim] = max(1, factor)
        if factors_for_level:
            inferred.append(factors_for_level)
    return inferred


def chunks_from_input(scale0: Any, label_dims: list[str]) -> dict[str, int] | None:
    chunksize = getattr(scale0.data, "chunksize", None)
    if not chunksize:
        return None
    return {
        dim: int(chunksize[idx])
        for idx, dim in enumerate(scale0.dims)
        if dim in label_dims
    }


def chunks_by_scale_from_input(multiscales: Any, label_dims: list[str]) -> dict[str, int] | None:
    images = getattr(multiscales, "images", None) or []
    if not images:
        return None

    # Use only scale-0 chunk sizes; to_multiscales expects one chunk spec shape,
    # not a per-pyramid-level sequence.
    first_image = images[0]
    chunksize = getattr(first_image.data, "chunksize", None)
    if not chunksize:
        return None

    dim_to_chunk = {
        dim: int(chunksize[idx])
        for idx, dim in enumerate(first_image.dims)
        if dim in label_dims
    }
    if any(dim not in dim_to_chunk for dim in label_dims):
        return None

    return dim_to_chunk


def resolve_cellpose_gpu_device() -> dict[str, Any]:
    """Require CUDA GPU and fail fast when GPU initialization is not possible."""
    try:
        import torch
    except ImportError as exc:
        raise RuntimeError("GPU_INIT_ERROR: torch is not installed in this environment") from exc

    try:
        cuda_available = bool(torch.cuda.is_available())
        cuda_device_count = int(torch.cuda.device_count()) if cuda_available else 0
    except Exception as exc:
        raise RuntimeError(f"GPU_INIT_ERROR: failed to query CUDA availability: {exc}") from exc

    cuda_visible_devices = os.getenv("CUDA_VISIBLE_DEVICES", "<unset>")
    if not cuda_available or cuda_device_count < 1:
        raise RuntimeError(
            "GPU_INIT_ERROR: CUDA GPU not available "
            f"(is_available={cuda_available}, device_count={cuda_device_count}, "
            f"CUDA_VISIBLE_DEVICES={cuda_visible_devices})"
        )

    try:
        _ = torch.cuda.get_device_name(0)
    except Exception as exc:
        raise RuntimeError(
            "GPU_INIT_ERROR: CUDA device discovery failed despite visible device(s): "
            f"{exc}"
        ) from exc

    if not core.use_gpu():
        raise RuntimeError(
            "GPU_INIT_ERROR: Cellpose reported no usable CUDA GPU despite CUDA availability"
        )

    return {"gpu": True}


def _select_timepoints(
    data: Any,
    dims: list[str],
    timepoint_indices: list[int] | None,
) -> tuple[list[int | None], int | None]:
    t_axis = dims.index("t") if "t" in dims else None
    if t_axis is None:
        if timepoint_indices is not None:
            raise ValueError("Input has no time dimension, but --timepoints was provided")
        return [None], t_axis

    n_timepoints = data.shape[t_axis]
    if timepoint_indices is None:
        selected = list(range(n_timepoints))
    else:
        selected = [int(t) for t in timepoint_indices]
        if not selected:
            raise ValueError("--timepoints was provided but no indices were specified")
        for t in selected:
            if t < 0 or t >= n_timepoints:
                raise ValueError(f"timepoint index {t} out of range for {n_timepoints} timepoint(s)")
    return selected, t_axis


def _select_channels(data: Any, dims: list[str], channel_indices: list[int] | None) -> list[int] | None:
    c_axis = dims.index("c") if "c" in dims else None
    if c_axis is None:
        if channel_indices is not None:
            # Allow [0] or 0 if it is the only channel provided
            if channel_indices == [0] or channel_indices == 0:
                return None
            raise ValueError(f"Input has no channel dimension, but --channels was provided: {channel_indices}")
        return None

    n_channels = data.shape[c_axis]
    if channel_indices is None:
        selected = list(range(n_channels))
    else:
        selected = [int(c) for c in channel_indices]
        if not selected:
            raise ValueError("--channels was provided but no channel indices were specified")
        for c in selected:
            if c < 0 or c >= n_channels:
                raise ValueError(f"channel index {c} out of range for {n_channels} channel(s)")
    return selected


def _parse_index_cli_values(values: list[str] | None, option_name: str) -> list[int] | None:
    if values is None:
        return None

    tokens: list[str] = []
    for raw in values:
        text = str(raw).strip()
        if not text:
            continue
        text = re.sub(r'^"+|"+$', '', text)
        text = re.sub(r"^'+|'+$", '', text)
        parts = [part for part in re.split(r"[;,\s]+", text) if part]
        tokens.extend(parts)

    if not tokens:
        return None

    try:
        return [int(token) for token in tokens]
    except ValueError as exc:
        raise ValueError(
            f"Invalid {option_name} value. Expected integer indices separated by semicolon, comma, or whitespace."
        ) from exc


def _parse_optional_float_cli(value: str) -> float | None:
    text = str(value).strip()
    if text.lower() == "none":
        return None
    return float(text)


def _parse_optional_int_cli(value: str) -> int | None:
    text = str(value).strip()
    if text.lower() == "none":
        return None
    return int(text)


def _parse_bool_cli(value: str) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no", "none"}:
        return False
    raise argparse.ArgumentTypeError(f"invalid boolean value: {value!r}")


def _normalize_bounds_arg(
    min_value: float | None,
    max_value: float | None,
    axis_name: str,
) -> tuple[float | None, float | None] | None:
    if min_value is None and max_value is None:
        return None

    start = float(min_value) if min_value is not None else None
    end = float(max_value) if max_value is not None else None
    if start is not None and end is not None and end <= start:
        raise ValueError(f"{axis_name}-bounds are invalid: min={start} must be < max={end}")
    return start, end


def _resolve_spatial_roi(
    data: Any,
    dims: list[str],
    scales: dict[str, float],
    translation: dict[str, float],
    x_min: float | None,
    x_max: float | None,
    y_min: float | None,
    y_max: float | None,
    z_min: float | None,
    z_max: float | None,
    halo: float = 0.0,
    use_physical_units: bool = False,
) -> dict[str, Any]:
    x_bounds = _normalize_bounds_arg(x_min, x_max, "x")
    y_bounds = _normalize_bounds_arg(y_min, y_max, "y")
    z_bounds = _normalize_bounds_arg(z_min, z_max, "z")

    if "x" not in dims or "y" not in dims:
        raise ValueError(f"Input dims must contain x and y axes, got {dims}")

    has_z_axis = "z" in dims
    if z_bounds is not None and not has_z_axis:
        raise ValueError("Input has no z dimension, but --z-min/--z-max was provided")

    if halo < 0.0:
        raise ValueError(f"--halo must be >= 0, got {halo}")

    resolved_bounds: dict[str, tuple[int, int]] = {}
    resolved_physical_bounds: dict[str, tuple[float, float]] = {}
    load_bounds: dict[str, tuple[int, int]] = {}
    load_physical_bounds: dict[str, tuple[float, float]] = {}
    has_user_roi = x_bounds is not None or y_bounds is not None or z_bounds is not None
    for dim in ("x", "y", "z"):
        if dim not in dims:
            continue

        axis = dims.index(dim)
        axis_size = int(data.shape[axis])
        scale = float(scales.get(dim, 1.0))
        if scale <= 0.0:
            raise ValueError(f"{dim} scale must be > 0, got {scale}")
        offset = float(translation.get(dim, 0.0))

        if dim == "x":
            interval = x_bounds
        elif dim == "y":
            interval = y_bounds
        else:
            interval = z_bounds

        if interval is None:
            start, end = 0, axis_size
            start_phys = offset
            end_phys = offset + (axis_size * scale)
        else:
            start_input, end_input = interval
            if use_physical_units:
                if start_input is None:
                    start = 0
                    start_phys = offset
                else:
                    start_phys = start_input
                    relative_start = (start_phys - offset) / scale
                    start = int(math.floor(relative_start + 1e-9))

                if end_input is None:
                    end = axis_size
                    end_phys = offset + (axis_size * scale)
                else:
                    end_phys = end_input
                    relative_end = (end_phys - offset) / scale
                    end = int(math.ceil(relative_end - 1e-9))

                if (start_input is not None and start < 0) or (end_input is not None and end > axis_size):
                    raise ValueError(
                        f"{dim}-bounds [{start_phys}, {end_phys}) are outside the image extent "
                        f"[{offset}, {offset + (axis_size * scale)}) at the selected resolution level"
                    )
            else:
                if start_input is None:
                    start = 0
                else:
                    start = int(math.floor(start_input + 1e-9))

                if end_input is None:
                    end = axis_size
                else:
                    end = int(math.ceil(end_input - 1e-9))

                if (start_input is not None and start < 0) or (end_input is not None and end > axis_size):
                    raise ValueError(
                        f"{dim}-bounds [{start_input}, {end_input}) in voxel units are outside the image extent "
                        f"[0, {axis_size}) at the selected resolution level"
                    )

                start_phys = offset + (start * scale)
                end_phys = offset + (end * scale)

            if end <= start:
                raise ValueError(
                    f"{dim}-bounds produce an empty ROI after conversion: start={start}, end={end}"
                )

        resolved_bounds[dim] = (start, end)
        resolved_physical_bounds[dim] = (start_phys, end_phys)

        if interval is None or halo == 0.0:
            load_start, load_end = start, end
        else:
            if use_physical_units:
                halo_voxels = int(math.ceil((halo / scale) - 1e-9))
            else:
                halo_voxels = int(math.ceil(halo - 1e-9))
            halo_voxels = max(0, halo_voxels)
            load_start = max(0, start - halo_voxels)
            load_end = min(axis_size, end + halo_voxels)

        load_bounds[dim] = (load_start, load_end)
        load_physical_bounds[dim] = (
            offset + (load_start * scale),
            offset + (load_end * scale),
        )

    load_slices = {dim: slice(start, end) for dim, (start, end) in load_bounds.items()}
    return {
        "has_user_roi": has_user_roi,
        "base_slices": {dim: slice(start, end) for dim, (start, end) in resolved_bounds.items()},
        "base_bounds": resolved_bounds,
        "base_physical_bounds": resolved_physical_bounds,
        "load_slices": load_slices,
        "load_bounds": load_bounds,
        "load_physical_bounds": load_physical_bounds,
    }


def _match_mask_shape(mask: np.ndarray, expected_shape: tuple[int, ...]) -> np.ndarray:
    """Match mask shape exactly to expected shape using edge-crop / zero-pad.

    The input image origin is preserved by applying all adjustments at the
    high-end (max index) of each axis.
    """
    if mask.ndim != len(expected_shape):
        raise ValueError(
            f"Cellpose mask rank mismatch: mask.ndim={mask.ndim}, expected ndim={len(expected_shape)}"
        )

    adjusted = mask
    for axis, target in enumerate(expected_shape):
        current = adjusted.shape[axis]
        if current > target:
            slicer = [slice(None)] * adjusted.ndim
            slicer[axis] = slice(0, target)
            adjusted = adjusted[tuple(slicer)]
        elif current < target:
            pad_width = [(0, 0)] * adjusted.ndim
            pad_width[axis] = (0, target - current)
            adjusted = np.pad(adjusted, pad_width=pad_width, mode="constant", constant_values=0)

    return np.asarray(adjusted)


def run_cellpose(
    data: Any,
    dims: list[str],
    scales: dict[str, float],
    translation: dict[str, float],
    cellpose_model_path: str = "",
    channel_indices: list[int] | None = None,
    timepoint_indices: list[int] | None = None,
    x_min: float | None = None,
    x_max: float | None = None,
    y_min: float | None = None,
    y_max: float | None = None,
    z_min: float | None = None,
    z_max: float | None = None,
    halo: float = 0.0,
    use_physical_units: bool = False,
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    do_3D: bool = False,
    flow3D_smooth: float = 0.0,
    stitch_threshold: float = 0.1,
    niter: int | None = None,
    min_size: int = 15,
) -> dict[str, Any]:
    model_kwargs = resolve_cellpose_gpu_device()
    if cellpose_model_path:
        model = models.CellposeModel(pretrained_model=cellpose_model_path, **model_kwargs)
    else:
        model = models.CellposeModel(**model_kwargs)

    selected_timepoints, t_axis = _select_timepoints(data, dims, timepoint_indices)
    selected_channels = _select_channels(data, dims, channel_indices)
    roi = _resolve_spatial_roi(
        data=data,
        dims=dims,
        scales=scales,
        translation=translation,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        z_min=z_min,
        z_max=z_max,
        halo=halo,
        use_physical_units=use_physical_units,
    )

    scale_x = float(scales.get("x", 1.0))
    scale_z = float(scales.get("z", 1.0))
    anisotropy = scale_z / scale_x if scale_x != 0.0 else 1.0

    label_mask_by_timepoint: list[np.ndarray] = []
    label_mask_dims_no_t: list[str] | None = None

    for t in selected_timepoints:
        indexer = [slice(None)] * data.ndim
        for dim, axis_slice in roi["load_slices"].items():
            indexer[dims.index(dim)] = axis_slice
        if t_axis is not None:
            indexer[t_axis] = t

        input_data = np.asarray(data[tuple(indexer)])
        input_dims = [dim for dim in dims if dim != "t"]

        if "c" in input_dims:
            c_axis_current = input_dims.index("c")
            if selected_channels is None:
                raise ValueError("No channel selection available despite channel dimension in input")
            if len(selected_channels) == 1:
                input_data = np.take(input_data, selected_channels[0], axis=c_axis_current)
                input_dims.pop(c_axis_current)
            else:
                input_data = np.take(input_data, selected_channels, axis=c_axis_current)

        if "z" in input_dims:
            n_slices = input_data.shape[input_dims.index("z")]
        else:
            n_slices = 1

        has_channel_axis = "c" in input_dims
        if has_channel_axis and input_data.ndim not in (3, 4):
            raise ValueError(f"Expected 2D/3D data with channel axis for Cellpose, got shape {input_data.shape}")
        if not has_channel_axis and input_data.ndim not in (2, 3):
            raise ValueError(f"Expected 2D or 3D data for Cellpose, got shape {input_data.shape}")

        eval_kwargs: dict[str, Any] = {
            "flow_threshold": flow_threshold,
            "cellprob_threshold": cellprob_threshold,
            "min_size": min_size,
            "niter": niter,
        }
        if diameter is not None:
            eval_kwargs["diameter"] = diameter
        if has_channel_axis:
            eval_kwargs["channel_axis"] = input_dims.index("c")
        if n_slices > 1:
            eval_kwargs["z_axis"] = input_dims.index("z")
            if do_3D:
                eval_kwargs["do_3D"] = True
                eval_kwargs["anisotropy"] = anisotropy
                eval_kwargs["flow3D_smooth"] = flow3D_smooth
            else:
                eval_kwargs["do_3D"] = False
                eval_kwargs["stitch_threshold"] = stitch_threshold

        label_mask, _, _ = model.eval(input_data, **eval_kwargs)

        # make sure that the cellpose SAM output has the exact same number
        # of pixels as the input
        expected_mask_dims = [dim for dim in input_dims if dim != "c"]
        expected_mask_shape = tuple(int(input_data.shape[input_dims.index(dim)]) for dim in expected_mask_dims)
        if tuple(label_mask.shape) != expected_mask_shape:
            print(
                "Adjusting Cellpose mask shape from "
                f"{tuple(label_mask.shape)} to {expected_mask_shape}"
            )
            label_mask = _match_mask_shape(np.asarray(label_mask), expected_mask_shape)

        # Segment with halo-expanded ROI but save labels for the originally requested ROI.
        crop_slicer = [slice(None)] * label_mask.ndim
        for axis, dim in enumerate(expected_mask_dims):
            if dim not in roi["base_bounds"]:
                continue
            base_start, base_end = roi["base_bounds"][dim]
            load_start, _ = roi["load_bounds"][dim]
            crop_start = base_start - load_start
            crop_end = crop_start + (base_end - base_start)
            crop_slicer[axis] = slice(crop_start, crop_end)
        label_mask = np.asarray(label_mask)[tuple(crop_slicer)]

        if label_mask.ndim == 2:
            current_label_mask_dims = ["y", "x"]
        elif label_mask.ndim == 3:
            current_label_mask_dims = ["z", "y", "x"]
        else:
            raise ValueError(f"Expected 2D or 3D masks, got shape {label_mask.shape}")

        if label_mask_dims_no_t is None:
            label_mask_dims_no_t = current_label_mask_dims
        elif label_mask_dims_no_t != current_label_mask_dims:
            raise ValueError(
                f"Inconsistent mask dimensions across timepoints: {label_mask_dims_no_t} vs {current_label_mask_dims}"
            )

        label_mask_by_timepoint.append(label_mask)

    if t_axis is not None:
        label_mask = np.stack(label_mask_by_timepoint, axis=0)
        if label_mask_dims_no_t is None:
            raise ValueError("Could not infer output mask dimensions")
        label_dims = ["t"] + label_mask_dims_no_t
    else:
        label_mask = label_mask_by_timepoint[0]
        if label_mask_dims_no_t is None:
            raise ValueError("Could not infer output mask dimensions")
        label_dims = label_mask_dims_no_t

    return {
        "label_mask": label_mask,
        "label_dims": label_dims,
        "model_kwargs": model_kwargs,
        "selected_timepoints": selected_timepoints,
        "selected_channels": selected_channels,
        "roi_bounds": roi["base_bounds"],
        "roi_physical_bounds": roi["base_physical_bounds"],
        "load_bounds": roi["load_bounds"],
        "load_physical_bounds": roi["load_physical_bounds"],
        "has_user_roi": roi["has_user_roi"],
        "halo": float(halo),
    }


def cellpose_zarr(
    input_zarr: str,
    output_zarr: str,
    cellpose_model_path: str = "",
    channel_indices: list[int] | None = None,
    timepoint_indices: list[int] | None = None,
    x_min: float | None = None,
    x_max: float | None = None,
    y_min: float | None = None,
    y_max: float | None = None,
    z_min: float | None = None,
    z_max: float | None = None,
    halo: float = 0.0,
    use_physical_units: bool = False,
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    do_3D: bool = False,
    flow3D_smooth: float = 0.0,
    stitch_threshold: float = 0.1,
    niter: int | None = None,
    min_size: int = 15,
    resolution_level: int = 0,
    ome_zarr_version: str | None = None,
) -> None:

    io.logger_setup()

    #
    # Open input OME-Zarr
    #

    print(f"Opening {input_zarr}")
    multiscales = nz.from_ngff_zarr(input_zarr)
    if resolution_level < 0:
        raise ValueError(f"--resolution-level must be >= 0, got {resolution_level}")
    if resolution_level >= len(multiscales.images):
        raise ValueError(
            f"--resolution-level {resolution_level} is out of range for input with {len(multiscales.images)} level(s)"
        )
    ngff_image = multiscales.images[resolution_level]
    data = ngff_image.data
    dims = list(ngff_image.dims)
    input_scales = dict(ngff_image.scale)
    input_axes_units = dict(getattr(ngff_image, "axes_units", {}) or {})
    input_translation = dict(getattr(ngff_image, "translation", {}) or {})
    
    #
    # Run CellPose
    #

    print("Running CellPose")
    result = run_cellpose(
        data=data,
        dims=dims,
        scales=input_scales,
        translation=input_translation,
        cellpose_model_path=cellpose_model_path,
        channel_indices=channel_indices,
        timepoint_indices=timepoint_indices,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        z_min=z_min,
        z_max=z_max,
        halo=halo,
        use_physical_units=use_physical_units,
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        do_3D=do_3D,
        flow3D_smooth=flow3D_smooth,
        stitch_threshold=stitch_threshold,
        niter=niter,
        min_size=min_size,
    )

    print(f"Cellpose device settings: {result['model_kwargs']}")
    print(f"Processed resolution level: {resolution_level}")
    print(f"Processed timepoint(s): {result['selected_timepoints']}")
    print(f"Used channel(s): {result['selected_channels']}")
    print(f"Requested ROI bounds (voxel): {result['roi_bounds']}")
    print(f"Requested ROI bounds (physical): {result['roi_physical_bounds']}")
    print(f"Loaded ROI bounds with halo (voxel): {result['load_bounds']}")
    print(f"Loaded ROI bounds with halo (physical): {result['load_physical_bounds']}")
    print(f"Applied halo: {result['halo']}")
    print(f"ROI input units: {'physical' if use_physical_units else 'voxel'}")


    #
    # Save output label mask as OME-Zarr
    #

    print(f"Writing labels to {output_zarr}")
    label_mask: np.ndarray = result["label_mask"]
    label_dims: list[str] = result["label_dims"]
    roi_bounds: dict[str, tuple[int, int]] = result["roi_bounds"]
    roi_physical_bounds: dict[str, tuple[float, float]] = result["roi_physical_bounds"]


    label_scales = {dim: float(input_scales.get(dim, 1.0)) for dim in label_dims}
    label_translation = {dim: float(input_translation.get(dim, 0.0)) for dim in label_dims}
    for dim in ("x", "y", "z"):
        if dim in label_dims and dim in roi_physical_bounds:
            label_translation[dim] = float(roi_physical_bounds[dim][0])
    label_axes_units = {dim: input_axes_units[dim] for dim in label_dims if dim in input_axes_units}

    label_image = nz.to_ngff_image(
        data=label_mask,
        name="nuclei_labels",
        dims=label_dims,
        scale=label_scales,
        translation=label_translation,
        axes_units=label_axes_units if label_axes_units else None,
    )

    input_scale_factors = scale_factors_from_input(multiscales, label_dims, resolution_level)
    chunks_spec: Any = chunks_from_input(ngff_image, label_dims)
    if chunks_spec is None:
        chunks_spec = None

    # Use nearest-neighbor pyramid generation for label data to avoid averaging label IDs.
    to_multiscales_kwargs: dict[str, Any] = {
        "scale_factors": input_scale_factors,
        "method": Methods.DASK_IMAGE_NEAREST,
    }
    if chunks_spec is not None:
        to_multiscales_kwargs["chunks"] = chunks_spec

    label_multiscales = nz.to_multiscales(label_image, **to_multiscales_kwargs)
    
    if ome_zarr_version is None:
        version = detect_ome_zarr_version(input_zarr)
    else:
        version = normalize_ome_zarr_version(ome_zarr_version)
    
    print(f"OME-Zarr version: {version}")
    print(f"OME-Zarr scale factors: {input_scale_factors}")
    print(f"OME-Zarr chunking: {chunks_spec}")
    print(f"OME-Zarr downsampling method: {to_multiscales_kwargs['method']}")

    nz.to_ngff_zarr(output_zarr, label_multiscales, version=version)

    # Cellpose instance IDs must be an integer dtype.
    label_root = zarr.open_group(output_zarr, mode="a")
    image_label_metadata = {
            "source": {
            "image": "../"
            }
        }

    if version == "0.5":
        ome = dict(label_root.attrs.get("ome", {}))
        ome["image-label"] = image_label_metadata
        label_root.attrs["ome"] = ome
    else:
        label_root.attrs["image-label"] = image_label_metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", 
        "--input_zarr", 
        type=str, 
        required=True)
    parser.add_argument(
        "-o", 
        "--output_zarr", 
        type=str, 
        required=True)
    parser.add_argument(
        "-c",
        "--channels",
        "--segmentation_channels",
        dest="segmentation_channels",
        type=str,
        nargs="+",
        default=None,
        help="Channel indices to process. Accepts semicolon/comma/space-separated values. If omitted, all channels are used.",
    )
    parser.add_argument(
        "-m", 
        "--cellpose_model_path", 
        type=str, 
        default=""
    )
    parser.add_argument(
        "-t",
        "--timepoints",
        type=int,
        nargs="+",
        default=None,
        help="Timepoint indices to process. If omitted, all timepoints are used."
    )
    parser.add_argument(
        "--ome-zarr-version",
        type=str,
        choices=["0.4", "0.5"],
        default=None,
        help="Target OME-Zarr version for output. If omitted, detected from input metadata."
    )
    parser.add_argument("--x-min", type=float, default=None, help="Optional x-axis ROI lower bound (inclusive).")
    parser.add_argument("--x-max", type=float, default=None, help="Optional x-axis ROI upper bound (exclusive).")
    parser.add_argument("--y-min", type=float, default=None, help="Optional y-axis ROI lower bound (inclusive).")
    parser.add_argument("--y-max", type=float, default=None, help="Optional y-axis ROI upper bound (exclusive).")
    parser.add_argument("--z-min", type=float, default=None, help="Optional z-axis ROI lower bound (inclusive).")
    parser.add_argument("--z-max", type=float, default=None, help="Optional z-axis ROI upper bound (exclusive).")
    parser.add_argument(
        "--halo",
        type=float,
        default=0.0,
        help="Optional non-negative ROI halo margin in the same units as --x-*/--y-*/--z-* (voxel by default, physical with --use-physical-units). Ignored when no ROI bounds are provided.",
    )
    parser.add_argument(
        "--use-physical-units",
        action="store_true",
        help="Interpret --x-*/--y-*/--z-* bounds in physical units instead of voxel units.",
    )
    parser.add_argument(
        "--resolution-level",
        type=int,
        default=0,
        help="Zero-based OME-Zarr resolution layer to operate on. Default: 0.",
    )
    parser.add_argument(
        "--diameter",
        type=_parse_optional_float_cli,
        default=None,
        help="Expected cell diameter in pixels, or 'None' to auto-estimate.",
    )
    parser.add_argument(
        "--flow_threshold",
        type=float,
        default=0.4,
        help="Cellpose flow error threshold.",
    )
    parser.add_argument(
        "--cellprob_threshold",
        type=float,
        default=0.0,
        help="Cellpose cell probability threshold.",
    )
    parser.add_argument(
        "--do_3D",
        type=_parse_bool_cli,
        default=False,
        help="Force 3D segmentation. If False (default), 3D mode is auto-selected from image anisotropy.",
    )
    parser.add_argument(
        "--flow3D_smooth",
        type=float,
        default=0.0,
        help="Sigma for smoothing flows in 3D segmentation.",
    )
    parser.add_argument(
        "--stitch_threshold",
        type=float,
        default=0.1,
        help="IoU threshold used to stitch per-slice 2D masks into 3D labels when running non-3D segmentation on a z-stack.",
    )
    parser.add_argument(
        "--niter",
        type=_parse_optional_int_cli,
        default=None,
        help="Number of mask reconstruction iterations, or 'None' to use Cellpose's default.",
    )
    parser.add_argument(
        "--min_size",
        type=int,
        default=15,
        help="Minimum number of pixels for a valid mask.",
    )

    args = parser.parse_args()
    try:
        segmentation_channels = _parse_index_cli_values(
            args.segmentation_channels,
            "--segmentation_channels/--channels",
        )
    except ValueError as exc:
        parser.error(str(exc))

    cellpose_zarr(
        args.input_zarr,
        args.output_zarr,
        args.cellpose_model_path,
        segmentation_channels,
        args.timepoints,
        args.x_min,
        args.x_max,
        args.y_min,
        args.y_max,
        args.z_min,
        args.z_max,
        args.halo,
        args.use_physical_units,
        diameter=args.diameter,
        flow_threshold=args.flow_threshold,
        cellprob_threshold=args.cellprob_threshold,
        do_3D=args.do_3D,
        flow3D_smooth=args.flow3D_smooth,
        stitch_threshold=args.stitch_threshold,
        niter=args.niter,
        min_size=args.min_size,
        resolution_level=args.resolution_level,
        ome_zarr_version=args.ome_zarr_version,
    )
