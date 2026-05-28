"""Per-cell QC metric computation for spatial transcriptomics.

Ported from R/QC.R (spatialPerCellQC) in SpaceTrooper.
"""

import warnings
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import anndata as ad
except ImportError:
    ad = None


def _find_negative_probes(feature_names: np.ndarray,
                           patterns: List[str]) -> dict:
    """Find negative probe indices by pattern matching on feature names.

    Parameters
    ----------
    feature_names : np.ndarray
        Array of feature (gene) names.
    patterns : list of str
        Patterns to match (prefix match).

    Returns
    -------
    dict
        Mapping of pattern to list of matching indices.
    """
    idx_map = {}
    for pattern in patterns:
        matches = [i for i, name in enumerate(feature_names)
                   if name.startswith(pattern)]
        if matches:
            idx_map[pattern] = matches
    return idx_map


def _compute_per_cell_qc(counts, subsets: dict) -> pd.DataFrame:
    """Compute per-cell QC metrics from count matrix.

    Parameters
    ----------
    counts : sparse or dense matrix
        Count matrix (cells x features).
    subsets : dict
        Mapping of subset name to feature indices.

    Returns
    -------
    DataFrame
        DataFrame with sum, detected, and subset metrics.
    """
    if sparse.issparse(counts):
        counts_dense = counts.toarray()
    else:
        counts_dense = np.asarray(counts)

    n_cells = counts_dense.shape[0]
    result = pd.DataFrame(index=range(n_cells))

    # Total counts and detected features per cell
    result["sum"] = counts_dense.sum(axis=1)
    result["detected"] = (counts_dense > 0).sum(axis=1)

    # Subset metrics
    for subset_name, idx_list in subsets.items():
        subset_counts = counts_dense[:, idx_list]
        safe_name = subset_name.replace(" ", "_").replace("-", "_")
        result[f"subsets_{safe_name}_sum"] = subset_counts.sum(axis=1)
        result[f"subsets_{safe_name}_detected"] = (subset_counts > 0).sum(axis=1)

    return result


def _compute_border_distance_cosmx(adata, xwindim: float,
                                    ywindim: float) -> None:
    """Compute minimum distance of each cell to the FOV border.

    Parameters
    ----------
    adata : AnnData
        Spatial data with FOV positions in .uns and spatial coords in .obs.
    xwindim : float
        FOV width in x.
    ywindim : float
        FOV height in y.
    """
    fov_positions = adata.uns.get("fov_positions")
    if fov_positions is None:
        return

    cd = adata.obs
    if "fov" not in cd.columns:
        return

    # Merge FOV positions
    if isinstance(fov_positions, pd.DataFrame):
        fov_df = fov_positions
    else:
        fov_df = pd.DataFrame(fov_positions)

    # Find coordinate columns (prefer global coordinates)
    spcn = [c for c in cd.columns if "global_px" in c.lower() and
            c.startswith(("CenterX", "CenterY", "x_", "y_"))]
    if len(spcn) < 2:
        spcn = [c for c in cd.columns if c.startswith(("CenterX", "CenterY",
                                                         "x_centroid", "y_centroid",
                                                         "center_x", "center_y"))]
    if len(spcn) < 2:
        spcn = [c for c in cd.columns if "x_global_px" in c or "y_global_px" in c]

    if len(spcn) < 2:
        return

    fovpn = [c for c in fov_df.columns if "x_global_px" in c or "y_global_px" in c]

    if len(fovpn) < 2 or "fov" not in fov_df.columns:
        return

    # Merge by fov
    fov_subset = fov_df[["fov"] + fovpn].copy()
    merged = cd.merge(fov_subset, on="fov", how="left")

    x_col = spcn[0]  # Cell X coordinate (e.g., CenterX_global_px)
    y_col = spcn[1]  # Cell Y coordinate (e.g., CenterY_global_px)
    x_fov_col = fovpn[0]  # FOV origin X (e.g., x_global_px)
    y_fov_col = fovpn[1]  # FOV origin Y (e.g., y_global_px)

    cd["dist_border_x"] = np.minimum(
        merged[x_col].values - merged[x_fov_col].values,
        (merged[x_fov_col].values + xwindim) - merged[x_col].values
    )
    cd["dist_border_y"] = np.minimum(
        merged[y_col].values - merged[y_fov_col].values,
        (merged[y_fov_col].values + ywindim) - merged[y_col].values
    )
    cd["dist_border"] = np.minimum(cd["dist_border_x"], cd["dist_border_y"])
    adata.obs = cd


def spatial_per_cell_qc(adata, micron_conv_fact: float = 0.12,
                         rm_zeros: bool = True,
                         neg_prob_list: Optional[List[str]] = None) -> None:
    """Compute quality-control metrics for each cell.

    Adds QC metrics to adata.obs. Equivalent to R spatialPerCellQC.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data.
    micron_conv_fact : float
        Factor to convert pixels to microns. Default 0.12.
    rm_zeros : bool
        Remove cells with zero counts. Default True.
    neg_prob_list : list of str, optional
        Patterns for negative probes. Default includes patterns for
        CosMx, Xenium, and MERFISH.
    """
    if neg_prob_list is None:
        neg_prob_list = [
            "NegPrb", "Negative", "SystemControl", "Ms IgG1", "Rb IgG",
            "BLANK_", "NegControlProbe", "NegControlCodeword",
            "UnassignedCodeword", "Blank",
        ]

    feature_names = np.array(adata.var_names)
    subsets = _find_negative_probes(feature_names, neg_prob_list)

    # Compute per-cell QC
    counts = adata.X
    qc_df = _compute_per_cell_qc(counts, subsets)

    # Copy QC columns to obs
    for col in qc_df.columns:
        adata.obs[col] = qc_df[col].values

    # Compute control metrics
    subset_sum_cols = [c for c in qc_df.columns if c.startswith("subsets_") and c.endswith("_sum")]
    subset_det_cols = [c for c in qc_df.columns if c.startswith("subsets_") and c.endswith("_detected")]

    if subset_sum_cols:
        npc = qc_df[subset_sum_cols].sum(axis=1).values
        npd = qc_df[subset_det_cols].sum(axis=1).values if subset_det_cols else np.zeros(len(qc_df))
    else:
        npc = np.zeros(len(qc_df))
        npd = np.zeros(len(qc_df))

    adata.obs["control_sum"] = npc
    adata.obs["control_detected"] = npd
    adata.obs["target_sum"] = adata.obs["sum"].values - npc
    adata.obs["target_detected"] = adata.obs["detected"].values - npd

    # Add spatial coordinates to obs if not already there
    if hasattr(adata, 'obsm') and 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
        if coords.shape[1] >= 2:
            coord_names = adata.uns.get("coord_names", ["x", "y"])
            for i, name in enumerate(coord_names[:2]):
                if name not in adata.obs.columns:
                    adata.obs[name] = coords[:, i]

    # Compute total (use 'total' from obs or compute from sum)
    if "total" not in adata.obs.columns:
        adata.obs["total"] = adata.obs["sum"]

    # Control-total ratio
    total = adata.obs["total"].values.astype(float)
    total_safe = np.where(total == 0, np.nan, total)
    adata.obs["ctrl_total_ratio"] = npc / total_safe
    adata.obs["ctrl_total_ratio"] = adata.obs["ctrl_total_ratio"].fillna(0)

    eps = 0.0001
    adata.obs["log2Ctrl_total_ratio"] = np.log2(
        adata.obs["ctrl_total_ratio"].values + eps
    )

    technology = adata.uns.get("technology", "")

    # CosMx-specific handling
    if technology in ("Nanostring_CosMx", "Nanostring_CosMx_Protein"):
        # Rename Area.um2 to Area_um if present
        if "Area.um2" in adata.obs.columns:
            adata.obs["Area_um"] = adata.obs["Area.um2"]

        # Convert coordinates from pixels to microns
        coord_cols = [c for c in adata.obs.columns
                     if "px" in c.lower() and ("global" in c.lower() or "center" in c.lower())]
        for col in coord_cols[:2]:
            new_name = col.replace("px", "um")
            adata.obs[new_name] = adata.obs[col] * micron_conv_fact

        # Compute Area_um
        if "Area_um" not in adata.obs.columns and "Area" in adata.obs.columns:
            adata.obs["Area_um"] = adata.obs["Area"] * (micron_conv_fact ** 2)

        # Compute border distance
        fov_dim = adata.uns.get("fov_dim", {"xdim": 4256, "ydim": 4256})
        _compute_border_distance_cosmx(adata, fov_dim["xdim"], fov_dim["ydim"])

    # Xenium-specific handling
    if technology == "10X_Xenium":
        if "Area_um" not in adata.obs.columns:
            if "cell_area" in adata.obs.columns:
                warnings.warn("Missing Area_um in colData for Xenium data.\n"
                              "Computing from cell_area column.")
                adata.obs["Area_um"] = adata.obs["cell_area"]

    # Aspect ratio
    if "AspectRatio" in adata.obs.columns:
        adata.obs["log2AspectRatio"] = np.log2(
            adata.obs["AspectRatio"].values.astype(float)
        )
    else:
        warnings.warn("Missing aspect ratio in colData")

    # Signal density
    if "Area_um" in adata.obs.columns:
        area = adata.obs["Area_um"].values.astype(float)
        area_safe = np.where(area == 0, np.nan, area)
        adata.obs["SignalDensity"] = adata.obs["sum"].values / area_safe
    else:
        adata.obs["SignalDensity"] = adata.obs["sum"].values

    if technology == "Nanostring_CosMx_Protein":
        adata.obs["SignalDensity"] = adata.obs["total"]

    adata.obs["log2SignalDensity"] = np.log2(
        adata.obs["SignalDensity"].values.astype(float)
    )

    # Remove zero-count cells
    if rm_zeros:
        zero_mask = adata.obs["sum"] == 0
        n_zeros = zero_mask.sum()
        if n_zeros > 0:
            print(f"Removing {n_zeros} cells with 0 counts!")
            adata._inplace_subset_obs(~zero_mask)
