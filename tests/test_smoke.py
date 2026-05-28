"""Smoke tests — verify the package imports and basic functionality."""

import numpy as np
import pandas as pd


def test_import():
    """Package imports without error."""
    import spacetrooper
    assert hasattr(spacetrooper, "__version__")
    assert spacetrooper.__version__ == "0.1.0"


def test_import_all_functions():
    """All public functions are importable."""
    from spacetrooper import (
        spatial_per_cell_qc,
        compute_spatial_outlier,
        compute_threshold_flags,
        compute_outliers_qc_score,
        check_outliers,
        get_fences_outlier,
        compute_qc_score,
        compute_train_df,
        compute_lambda,
        train_model,
        get_model_formula,
        get_active_geometry_name,
        set_active_geometry,
        rename_geometry,
        read_polygons,
        read_polygons_cosmx,
        read_polygons_xenium,
        read_polygons_merfish,
        add_polygons_to_adata,
        compute_area_from_polygons,
        compute_aspect_ratio_from_polygons,
        compute_center_from_polygons,
    )


def test_spatial_per_cell_qc_runs(tiny_cosmx_adata):
    """spatial_per_cell_qc runs without error."""
    from spacetrooper.qc_metrics import spatial_per_cell_qc

    adata = tiny_cosmx_adata.copy()
    spatial_per_cell_qc(adata, rm_zeros=False)

    # Check that QC metrics were added
    assert "control_sum" in adata.obs.columns
    assert "target_sum" in adata.obs.columns
    assert "ctrl_total_ratio" in adata.obs.columns
    assert "log2Ctrl_total_ratio" in adata.obs.columns
    assert "SignalDensity" in adata.obs.columns
    assert "log2SignalDensity" in adata.obs.columns


def test_compute_threshold_flags_runs(tiny_cosmx_with_qc):
    """compute_threshold_flags runs without error."""
    from spacetrooper.outlier_detection import compute_threshold_flags

    adata = tiny_cosmx_with_qc.copy()
    compute_threshold_flags(adata)

    assert "is_zero_counts" in adata.obs.columns
    assert "is_ctrl_tot_outlier" in adata.obs.columns
    assert "threshold_flags" in adata.obs.columns


def test_compute_spatial_outlier_runs(tiny_cosmx_with_qc):
    """compute_spatial_outlier runs without error."""
    from spacetrooper.outlier_detection import compute_spatial_outlier

    adata = tiny_cosmx_with_qc.copy()
    # Use scuttle method (MAD-based) since mc may fail on small/symmetric data
    compute_spatial_outlier(adata, compute_by="log2SignalDensity", method="scuttle")

    assert "log2SignalDensity_outlier_sc" in adata.obs.columns
    assert "log2SignalDensity_outlier_sc_thresholds" in adata.uns


def test_mc_computation():
    """Medcouple computation works on known data."""
    from spacetrooper.outlier_detection import _mc

    # Right-skewed data: mc should be positive
    rng = np.random.RandomState(42)
    x = rng.exponential(1, 1000)
    mc_val = _mc(x)
    assert mc_val > 0

    # Symmetric data: mc should be near 0
    x = rng.normal(0, 1, 1000)
    mc_val = _mc(x)
    assert abs(mc_val) < 0.2


def test_scuttle_outlier():
    """MAD-based outlier detection works."""
    from spacetrooper.outlier_detection import _is_outlier_scuttle

    rng = np.random.RandomState(42)
    x = rng.normal(0, 1, 1000)
    # Add some clear outliers
    x = np.append(x, [100, -100])

    outlier, (lower, upper) = _is_outlier_scuttle(x, type_="both")
    assert outlier.sum() >= 2  # At least the two clear outliers
    assert lower < upper


def test_parity_metrics_available():
    """Parity metrics module is accessible."""
    import sys
    sys.path.insert(0, "D:/test/omicverse-rebuildr-main/omicverse-rebuildr-main")
    from engine.parity_metrics import (
        parity_deterministic,
        parity_classification,
        is_pass,
    )
