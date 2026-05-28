"""Pytest fixtures for SpaceTrooper tests."""

import sys
from pathlib import Path

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

# Add package to path
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def tiny_cosmx_adata():
    """Create a minimal AnnData mimicking CosMx data for testing."""
    np.random.seed(42)
    n_cells = 100
    n_genes = 50

    # Generate count matrix
    counts = sparse.random(n_cells, n_genes, density=0.3, format="csr",
                          random_state=42)
    counts.data = np.round(counts.data * 100).astype(int)

    # Generate metadata
    obs = pd.DataFrame(index=[f"f1_c{i}" for i in range(n_cells)])
    obs["cell_id"] = obs.index
    obs["fov"] = 1
    obs["cellID"] = range(n_cells)
    obs["total"] = np.array(counts.sum(axis=1)).flatten()
    obs["sum"] = obs["total"]
    obs["detected"] = np.array((counts > 0).sum(axis=1)).flatten()
    obs["Area"] = np.random.uniform(10, 100, n_cells)
    obs["Area_um"] = obs["Area"] * 0.12**2
    obs["AspectRatio"] = np.random.uniform(0.5, 2.0, n_cells)

    # Generate spatial coordinates
    coords = np.column_stack([
        np.random.uniform(0, 4256, n_cells),
        np.random.uniform(0, 4256, n_cells),
    ])

    # Create gene names with some negative probes
    gene_names = [f"Gene_{i}" for i in range(n_genes - 5)]
    gene_names += ["NegPrb_1", "NegPrb_2", "Negative_1", "SystemControl_1", "BLANK_1"]

    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obsm["spatial"] = coords
    adata.uns["technology"] = "Nanostring_CosMx"
    adata.uns["sample_id"] = "test_sample"
    adata.uns["fov_dim"] = {"xdim": 4256, "ydim": 4256}
    adata.uns["fov_positions"] = pd.DataFrame({
        "fov": [1],
        "x_global_px": [0],
        "y_global_px": [0],
    })

    return adata


@pytest.fixture
def tiny_cosmx_with_qc(tiny_cosmx_adata):
    """CosMx AnnData with QC metrics already computed."""
    from spacetrooper.qc_metrics import spatial_per_cell_qc
    adata = tiny_cosmx_adata.copy()
    # Don't remove zeros to keep consistent cell count
    spatial_per_cell_qc(adata, rm_zeros=False)
    return adata


@pytest.fixture
def r_reference_dir():
    """Path to R reference source."""
    return Path(__file__).parent.parent / "SpaceTrooper-ref"
