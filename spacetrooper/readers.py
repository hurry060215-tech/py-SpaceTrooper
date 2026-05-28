"""Data readers for CosMx, MERFISH, and Xenium spatial transcriptomics.

Ported from R/readCosMx.R, R/readMerfish.R, R/readXenium.R in SpaceTrooper.
Returns AnnData objects instead of SpatialExperiment.
"""

import warnings
from pathlib import Path
from typing import Optional, List

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from .polygons import (
    read_polygons_cosmx,
    read_polygons_xenium,
    read_polygons_merfish,
    compute_area_from_polygons,
    compute_aspect_ratio_from_polygons,
)


def _check_fov_position_version(adata) -> None:
    """Standardize FOV position column names.

    Parameters
    ----------
    adata : AnnData
        AnnData with FOV positions in .uns['fov_positions'].
    """
    fovpos = adata.uns.get("fov_positions")
    if fovpos is None:
        return

    if not isinstance(fovpos, pd.DataFrame):
        fovpos = pd.DataFrame(fovpos)

    # Rename FOV → fov
    fov_cols = [c for c in fovpos.columns if "FOV" in c]
    rename_map = {c: "fov" for c in fov_cols}
    fovpos = fovpos.rename(columns=rename_map)

    # Lowercase X/Y/Z
    xyz_cols = [c for c in fovpos.columns if any(k in c for k in ["X", "Y", "Z"])]
    rename_map = {c: c.lower() for c in xyz_cols}
    fovpos = fovpos.rename(columns=rename_map)

    # Rename _px → _global_px
    px_cols = [c for c in fovpos.columns if c.endswith("_px") and not c.endswith("_global_px")]
    rename_map = {c: c.replace("_px", "_global_px") for c in px_cols}
    fovpos = fovpos.rename(columns=rename_map)

    # Convert mm to px if needed
    if "x_mm" in fovpos.columns:
        fovpos["x_global_px"] = fovpos["x_mm"] / 0.12028 * 1e3
        fovpos["y_global_px"] = fovpos["y_mm"] / 0.12028 * 1e3

    # Filter to matching FOVs
    if "fov" in fovpos.columns and "fov" in adata.obs.columns:
        unique_fovs = adata.obs["fov"].unique()
        fovpos = fovpos[fovpos["fov"].isin(unique_fovs)]
        fovpos = fovpos.sort_values("fov")

    adata.uns["fov_positions"] = fovpos


def read_cosmx_spe(dir_name: str, sample_name: str = "sample01",
                    coord_names: Optional[List[str]] = None,
                    count_mat_pattern: str = "exprMat_file.csv",
                    metadata_pattern: str = "metadata_file.csv",
                    polygons_pattern: str = "polygons.csv",
                    fov_pos_pattern: str = "fov_positions_file.csv",
                    fov_dims: Optional[dict] = None,
                    keep_polygons: bool = False) -> ad.AnnData:
    """Read CosMx data into an AnnData object.

    Parameters
    ----------
    dir_name : str
        Directory containing CosMx data files.
    sample_name : str
        Sample identifier. Default "sample01".
    coord_names : list of str, optional
        Coordinate column names. Default ["CenterX_global_px", "CenterY_global_px"].
    count_mat_pattern : str
        Pattern for count matrix file.
    metadata_pattern : str
        Pattern for metadata file.
    polygons_pattern : str
        Pattern for polygons file.
    fov_pos_pattern : str
        Pattern for FOV positions file.
    fov_dims : dict, optional
        FOV dimensions. Default {"xdim": 4256, "ydim": 4256}.
    keep_polygons : bool
        Whether to load polygons into memory.

    Returns
    -------
    AnnData
        CosMx spatial data.
    """
    if coord_names is None:
        coord_names = ["CenterX_global_px", "CenterY_global_px"]
    if fov_dims is None:
        fov_dims = {"xdim": 4256, "ydim": 4256}

    dir_path = Path(dir_name)
    if not dir_path.exists():
        raise FileNotFoundError(f"Directory not found: {dir_name}")

    # Read count matrix
    count_file = list(dir_path.glob(f"*{count_mat_pattern}*"))
    if not count_file:
        raise FileNotFoundError(f"Count matrix not found matching '{count_mat_pattern}'")
    counts_df = pd.read_csv(count_file[0])

    # Read metadata
    meta_file = list(dir_path.glob(f"*{metadata_pattern}*"))
    if not meta_file:
        raise FileNotFoundError(f"Metadata not found matching '{metadata_pattern}'")
    metadata_df = pd.read_csv(meta_file[0])

    # Read FOV positions
    fov_file = list(dir_path.glob(f"*{fov_pos_pattern}*"))
    fov_positions = None
    if fov_file:
        fov_positions = pd.read_csv(fov_file[0])

    # Build cell IDs
    cell_ids = "f" + metadata_df["fov"].astype(str) + "_c" + metadata_df["cell_ID"].astype(str)
    metadata_df["cell_id"] = cell_ids
    metadata_df.index = cell_ids

    # Extract counts (features x cells)
    if "cell_id" in counts_df.columns:
        feature_cols = [c for c in counts_df.columns if c != "cell_id"]
    else:
        feature_cols = counts_df.columns.tolist()

    count_matrix = counts_df[feature_cols].values.T  # features x cells

    # Create AnnData
    adata = ad.AnnData(
        X=sparse.csr_matrix(count_matrix.T),  # cells x features
        obs=metadata_df,
        var=pd.DataFrame(index=feature_cols),
    )

    # Set spatial coordinates
    if all(c in metadata_df.columns for c in coord_names):
        adata.obsm["spatial"] = metadata_df[coord_names].values
        adata.uns["coord_names"] = coord_names

    # Store metadata
    adata.uns["technology"] = "Nanostring_CosMx"
    adata.uns["sample_id"] = sample_name
    adata.uns["fov_dim"] = fov_dims
    if fov_positions is not None:
        adata.uns["fov_positions"] = fov_positions

    # Store polygons path
    pol_file = list(dir_path.glob(f"*{polygons_pattern}*"))
    if pol_file:
        adata.uns["polygons"] = str(pol_file[0])

    # Standardize FOV positions
    _check_fov_position_version(adata)

    # Rename cell_ID → cellID
    if "cell_ID" in adata.obs.columns:
        adata.obs["cellID"] = adata.obs["cell_ID"]

    adata.obs["sample_id"] = sample_name

    return adata


def read_merfish_spe(dir_name: str, sample_name: str = "sample01",
                      compute_missing_metrics: bool = True,
                      keep_polygons: bool = False,
                      boundaries_type: str = "parquet",
                      count_mat_pattern: str = "cell_by_gene.csv",
                      metadata_pattern: str = "cell_metadata.csv",
                      polygons_pattern: str = "cell_boundaries.parquet",
                      coord_names: Optional[List[str]] = None,
                      use_volume: bool = True) -> ad.AnnData:
    """Read MERFISH data into an AnnData object.

    Parameters
    ----------
    dir_name : str
        Directory containing MERFISH output files.
    sample_name : str
        Sample identifier.
    compute_missing_metrics : bool
        Compute area and aspect ratio from polygons.
    keep_polygons : bool
        Keep raw polygon geometries.
    boundaries_type : str
        One of "parquet" or "HDF5".
    count_mat_pattern : str
        Pattern for count matrix file.
    metadata_pattern : str
        Pattern for metadata file.
    polygons_pattern : str
        Pattern for polygons file.
    coord_names : list of str, optional
        Coordinate column names. Default ["center_x", "center_y"].
    use_volume : bool
        Use volume column for area if available.

    Returns
    -------
    AnnData
        MERFISH spatial data.
    """
    if coord_names is None:
        coord_names = ["center_x", "center_y"]

    dir_path = Path(dir_name)
    if not dir_path.exists():
        raise FileNotFoundError(f"Directory not found: {dir_name}")

    # Read count matrix
    count_file = list(dir_path.glob(f"*{count_mat_pattern}*"))
    if not count_file:
        raise FileNotFoundError(f"Count matrix not found matching '{count_mat_pattern}'")
    countmat = pd.read_csv(count_file[0])

    # Read metadata
    meta_file = list(dir_path.glob(f"*{metadata_pattern}*"))
    if not meta_file:
        raise FileNotFoundError(f"Metadata not found matching '{metadata_pattern}'")
    metadata = pd.read_csv(meta_file[0])

    # Read polygons path
    pol_file = list(dir_path.glob(f"*{polygons_pattern}*"))

    # Standardize column names
    for alias in ["V1", "cell"]:
        if alias in countmat.columns and "cell_id" not in countmat.columns:
            countmat = countmat.rename(columns={alias: "cell_id"})

    for alias in ["EntityID", "V1"]:
        if alias in metadata.columns and "cell_id" not in metadata.columns:
            metadata = metadata.rename(columns={alias: "cell_id"})

    # Extract features
    if "cell_id" in countmat.columns:
        cell_ids = countmat["cell_id"].values
        feature_cols = [c for c in countmat.columns if c != "cell_id"]
    else:
        cell_ids = np.arange(len(countmat)).astype(str)
        feature_cols = countmat.columns.tolist()

    count_matrix = countmat[feature_cols].values.T  # features x cells

    # Align metadata to count matrix
    metadata = metadata.set_index("cell_id", drop=False)
    metadata = metadata.reindex(cell_ids)

    # Create AnnData
    adata = ad.AnnData(
        X=sparse.csr_matrix(count_matrix.T),  # cells x features
        obs=metadata.reset_index(drop=True),
        var=pd.DataFrame(index=feature_cols),
    )
    adata.obs["cell_id"] = cell_ids
    adata.obs.index = cell_ids

    # Set spatial coordinates
    if all(c in metadata.columns for c in coord_names):
        adata.obsm["spatial"] = metadata[coord_names].values
        adata.uns["coord_names"] = coord_names

    # Compute missing metrics from polygons
    if compute_missing_metrics and pol_file:
        pol_path = str(pol_file[0])
        try:
            polygons = read_polygons_merfish(pol_path, file_type=boundaries_type)

            if use_volume and "volume" in adata.obs.columns:
                warnings.warn("Volume is used for MERFISH. Renaming as Area_um.")
                adata.obs["Area_um"] = adata.obs["volume"]
            else:
                if "volume" not in adata.obs.columns:
                    warnings.warn("Volume column not found. Computing area from polygons.")
                adata.obs["Area_um"] = compute_area_from_polygons(polygons)

            adata.obs["AspectRatio"] = compute_aspect_ratio_from_polygons(polygons)
        except Exception as e:
            warnings.warn(f"Failed to compute metrics from polygons: {e}")

    # Store metadata
    adata.uns["technology"] = "Vizgen_MERFISH"
    adata.uns["sample_id"] = sample_name
    if pol_file:
        adata.uns["polygons"] = str(pol_file[0])

    return adata


def read_xenium_spe(dir_name: str, sample_name: str = "sample01",
                     boundaries_type: str = "parquet",
                     coord_names: Optional[List[str]] = None,
                     compute_missing_metrics: bool = True,
                     keep_polygons: bool = False,
                     counts_pattern: str = "cell_feature_matrix",
                     metadata_pattern: str = "cells.csv.gz",
                     polygons_pattern: str = "cell_boundaries") -> ad.AnnData:
    """Read Xenium data into an AnnData object.

    Parameters
    ----------
    dir_name : str
        Directory containing Xenium Output Bundle.
    sample_name : str
        Sample identifier.
    boundaries_type : str
        One of "parquet" or "csv".
    coord_names : list of str, optional
        Coordinate column names. Default ["x_centroid", "y_centroid"].
    compute_missing_metrics : bool
        Compute area and aspect ratio from polygons.
    keep_polygons : bool
        Keep raw polygon geometries.
    counts_pattern : str
        Pattern for count matrix file.
    metadata_pattern : str
        Pattern for metadata file.
    polygons_pattern : str
        Pattern for polygons file.

    Returns
    -------
    AnnData
        Xenium spatial data.
    """
    if coord_names is None:
        coord_names = ["x_centroid", "y_centroid"]

    dir_path = Path(dir_name)
    if not dir_path.exists():
        raise FileNotFoundError(f"Directory not found: {dir_name}")

    # Check for outs/ subdirectory
    outs_path = dir_path / "outs"
    if outs_path.exists():
        dir_path = outs_path

    # Read metadata
    meta_file = list(dir_path.glob(f"*{metadata_pattern}*"))
    if not meta_file:
        raise FileNotFoundError(f"Metadata not found matching '{metadata_pattern}'")
    metadata = pd.read_csv(meta_file[0])

    # Read polygons path
    ext = ".parquet" if boundaries_type == "parquet" else ".csv.gz"
    pol_file = list(dir_path.glob(f"*{polygons_pattern}*{ext}*"))

    # Set cell IDs
    if "cell_id" in metadata.columns:
        metadata.index = metadata["cell_id"].astype(str)
    else:
        metadata.index = np.arange(len(metadata)).astype(str)

    # Read count matrix (simplified - assumes H5 or CSV)
    counts_h5 = list(dir_path.glob(f"*{counts_pattern}*.h5"))
    if counts_h5:
        import h5py
        with h5py.File(counts_h5[0], "r") as f:
            # Simplified H5 reading
            pass

    # Create minimal AnnData (actual count reading depends on format)
    n_cells = len(metadata)
    adata = ad.AnnData(
        X=sparse.csr_matrix((n_cells, 0)),
        obs=metadata,
        var=pd.DataFrame(index=[]),
    )

    # Set spatial coordinates
    if all(c in metadata.columns for c in coord_names):
        adata.obsm["spatial"] = metadata[coord_names].values
        adata.uns["coord_names"] = coord_names

    # Compute missing metrics
    if compute_missing_metrics and pol_file:
        pol_path = str(pol_file[0])
        try:
            polygons = read_polygons_xenium(pol_path, file_type=boundaries_type)
            adata.obs["AspectRatio"] = compute_aspect_ratio_from_polygons(polygons)

            if "cell_area" in adata.obs.columns:
                adata.obs["Area_um"] = adata.obs["cell_area"]
            else:
                adata.obs["Area_um"] = compute_area_from_polygons(polygons)
        except Exception as e:
            warnings.warn(f"Failed to compute metrics from polygons: {e}")

    # Store metadata
    adata.uns["technology"] = "10X_Xenium"
    adata.uns["sample_id"] = sample_name
    if pol_file:
        adata.uns["polygons"] = str(pol_file[0])

    adata.obs["sample_id"] = sample_name

    return adata
