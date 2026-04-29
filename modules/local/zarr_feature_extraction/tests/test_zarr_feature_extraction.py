import numpy as np
import pytest

from zarr_feature_extraction import build_label_path, parse_properties, prepare_intensity_image


def test_build_label_path_default():
    assert build_label_path("nuclei") == "labels/nuclei"


def test_build_label_path_explicit():
    assert build_label_path("labels/membrane") == "labels/membrane"


def test_parse_properties():
    assert parse_properties("label, area,centroid") == ["label", "area", "centroid"]


def test_prepare_intensity_image_moves_channel():
    intensity = np.zeros((1, 4, 4), dtype=np.uint8)
    out = prepare_intensity_image(intensity, (4, 4))
    assert out.shape == (4, 4, 1)


def test_prepare_intensity_image_shape_mismatch():
    intensity = np.zeros((2, 4, 4), dtype=np.uint8)
    with pytest.raises(ValueError, match="incompatible"):
        prepare_intensity_image(intensity, (5, 5))
