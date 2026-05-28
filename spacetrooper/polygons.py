"""Polygon reading and manipulation functions.

Ported from R/polygons.R in SpaceTrooper.
Uses geopandas/shapely instead of R sf.
"""

import warnings
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid


def _create_polygons(spat_obj: pd.DataFrame, x: str, y: str,
                     polygon_id: str = "cell_id") -> gpd.GeoDataFrame:
    """Create a GeoDataFrame of polygons from vertex-level data.

    Groups by polygon_id, collects (x, y) vertices, and builds Shapely Polygons.

    Parameters
    ----------
    spat_obj : DataFrame
        Must contain columns x, y, and polygon_id.
    x : str
        Column name for x-coordinates.
    y : str
        Column name for y-coordinates.
    polygon_id : str
        Column name for polygon identifiers.

    Returns
    -------
    GeoDataFrame
        One row per polygon with geometry column.
    """
    geometries = []
    ids = []
    extra_cols = {}

    for pid, group in spat_obj.groupby(polygon_id):
        coords = list(zip(group[x].values, group[y].values))
        if len(coords) < 3:
            continue
        # Close the polygon if not already closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        try:
            poly = Polygon(coords)
            if not poly.is_valid:
                poly = make_valid(poly)
            geometries.append(poly)
            ids.append(pid)
        except Exception:
            continue

    result = gpd.GeoDataFrame({polygon_id: ids}, geometry=geometries)
    result = result.rename(columns={polygon_id: "cell_id"})
    return result


def _check_polygons_validity(gdf: gpd.GeoDataFrame,
                              keep_multi_pol: bool = True,
                              verbose: bool = False) -> gpd.GeoDataFrame:
    """Check polygon validity and optionally remove multipolygons.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame with polygon geometries.
    keep_multi_pol : bool
        If False, remove multipolygons. Default True.
    verbose : bool
        Print progress messages.

    Returns
    -------
    GeoDataFrame
        Validated GeoDataFrame.
    """
    # Fix invalid geometries
    invalid_mask = ~gdf.geometry.is_valid
    if invalid_mask.any():
        gdf.loc[invalid_mask, "geometry"] = gdf.loc[invalid_mask, "geometry"].apply(
            lambda g: make_valid(g).buffer(0)
        )

    # Track multipolygons
    is_multi = gdf.geometry.apply(lambda g: isinstance(g, MultiPolygon))
    gdf = gdf.copy()
    gdf["is_multi"] = is_multi
    gdf["multi_n"] = gdf.geometry.apply(
        lambda g: len(g.geoms) if isinstance(g, MultiPolygon) else 1
    )

    if verbose:
        print(f"Detected {is_multi.sum()} multipolygons.")

    if not keep_multi_pol:
        if verbose:
            print(f"Removing {is_multi.sum()} multipolygons.")
        gdf = gdf[~gdf["is_multi"]].copy()

    return gdf


def read_polygons(polygons_file: str, file_type: str = "csv",
                  x: str = "x_global_px", y: str = "y_global_px",
                  xloc: Optional[str] = "x_local_px",
                  yloc: Optional[str] = "y_local_px",
                  keep_multi_pol: bool = True,
                  verbose: bool = False) -> gpd.GeoDataFrame:
    """Read and validate polygons from a file.

    Parameters
    ----------
    polygons_file : str
        Path to the polygon file.
    file_type : str
        One of "csv", "parquet". Default "csv".
    x : str
        Column name for x-coordinates. Default "x_global_px".
    y : str
        Column name for y-coordinates. Default "y_global_px".
    xloc : str, optional
        Column name for local x-coordinates.
    yloc : str, optional
        Column name for local y-coordinates.
    keep_multi_pol : bool
        Whether to keep multipolygons. Default True.
    verbose : bool
        Print progress messages.

    Returns
    -------
    GeoDataFrame
        Validated polygon geometries.
    """
    path = Path(polygons_file)
    if not path.exists():
        raise FileNotFoundError(f"Polygon file not found: {polygons_file}")

    if file_type == "csv":
        spat_obj = pd.read_csv(polygons_file)
    elif file_type == "parquet":
        spat_obj = pd.read_parquet(polygons_file)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    # Build cell_id if missing
    if "cell_id" not in spat_obj.columns:
        spat_obj["cell_id"] = "f" + spat_obj["fov"].astype(str) + "_c" + spat_obj["cellID"].astype(str)
    spat_obj["cell_id"] = spat_obj["cell_id"].astype(str)

    # Remove polygons with fewer than 4 points
    counts = spat_obj.groupby("cell_id").size()
    small_polys = counts[counts < 4]
    if len(small_polys) > 0:
        warnings.warn(f"Removing {len(small_polys)} polygons with less than 4 points.")
        valid_ids = counts[counts >= 4].index
        spat_obj = spat_obj[spat_obj["cell_id"].isin(valid_ids)]

    # Create global polygons
    polygons = _create_polygons(spat_obj, x=x, y=y, polygon_id="cell_id")
    polygons = rename_geometry(polygons, "geometry", "global")

    # Create local polygons if columns exist
    if xloc and yloc and xloc in spat_obj.columns and yloc in spat_obj.columns:
        polygons_loc = _create_polygons(spat_obj, x=xloc, y=yloc, polygon_id="cell_id")
        polygons["local"] = polygons_loc["geometry"]

    if verbose:
        print(f"Polygons detected: {len(polygons)}")

    polygons = _check_polygons_validity(polygons, keep_multi_pol=keep_multi_pol,
                                         verbose=verbose)
    polygons.index = polygons["cell_id"]

    if verbose:
        print(f"Polygons after validity: {len(polygons)}")

    return polygons


def read_polygons_cosmx(polygons_file: str, file_type: str = "csv",
                        x: str = "x_global_px", y: str = "y_global_px",
                        xloc: str = "x_local_px", yloc: str = "y_local_px",
                        keep_multi_pol: bool = True,
                        verbose: bool = False) -> gpd.GeoDataFrame:
    """Read polygon data specific to CosMx technology.

    Parameters
    ----------
    polygons_file : str
        Path to the polygon file.
    file_type : str
        One of "csv" or "parquet".
    x, y : str
        Column names for global coordinates.
    xloc, yloc : str
        Column names for local coordinates.
    keep_multi_pol : bool
        Whether to keep multipolygons.
    verbose : bool
        Print progress messages.

    Returns
    -------
    GeoDataFrame
        CosMx polygon data with global and local geometries.
    """
    polygons = read_polygons(polygons_file, file_type=file_type, x=x, y=y,
                             xloc=xloc, yloc=yloc, keep_multi_pol=keep_multi_pol,
                             verbose=verbose)
    # Reorder columns: mandatory first
    mandatory = ["cell_id", "global", "is_multi", "multi_n"]
    extra = [c for c in polygons.columns if c not in mandatory]
    polygons = polygons[mandatory + extra]
    return polygons


def read_polygons_xenium(polygons_file: str, file_type: str = "parquet",
                         x: str = "vertex_x", y: str = "vertex_y",
                         keep_multi_pol: bool = True,
                         verbose: bool = False) -> gpd.GeoDataFrame:
    """Read polygon data specific to Xenium technology.

    Parameters
    ----------
    polygons_file : str
        Path to the polygon file.
    file_type : str
        One of "parquet" or "csv".
    x, y : str
        Column names for vertex coordinates.
    keep_multi_pol : bool
        Whether to keep multipolygons.
    verbose : bool
        Print progress messages.

    Returns
    -------
    GeoDataFrame
        Xenium polygon data.
    """
    polygons = read_polygons(polygons_file, file_type=file_type, x=x, y=y,
                             xloc=None, yloc=None, keep_multi_pol=keep_multi_pol,
                             verbose=verbose)
    mandatory = ["cell_id", "global", "is_multi", "multi_n"]
    extra = [c for c in polygons.columns if c not in mandatory]
    polygons = polygons[mandatory + extra]
    return polygons


def read_polygons_merfish(polygons_source: str, file_type: str = "parquet",
                          z_lev: int = 3, z_column: str = "ZIndex",
                          keep_multi_pol: bool = True,
                          verbose: bool = False) -> gpd.GeoDataFrame:
    """Read polygon data specific to MERFISH technology.

    Parameters
    ----------
    polygons_source : str
        Path to parquet file or HDF5 folder.
    file_type : str
        One of "parquet" or "HDF5".
    z_lev : int
        Z-level to filter. Default 3.
    z_column : str
        Column name for Z index. Default "ZIndex".
    keep_multi_pol : bool
        Whether to keep multipolygons.
    verbose : bool
        Print progress messages.

    Returns
    -------
    GeoDataFrame
        MERFISH polygon data.
    """
    if file_type == "HDF5":
        raise NotImplementedError("HDF5 polygon reading not yet implemented for Python port")

    polygons = pd.read_parquet(polygons_source)
    polygons = polygons[polygons[z_column] == z_lev].copy()
    polygons["cell_id"] = polygons["EntityID"].astype(str)

    gdf = _create_polygons(polygons, x="vertex_x", y="vertex_y", polygon_id="cell_id")
    gdf = rename_geometry(gdf, "geometry", "global")
    gdf = _check_polygons_validity(gdf, keep_multi_pol=keep_multi_pol, verbose=verbose)

    mandatory = ["cell_id", "global", "is_multi", "multi_n"]
    extra = [c for c in gdf.columns if c not in mandatory]
    gdf = gdf[mandatory + extra]
    return gdf


def compute_area_from_polygons(gdf: gpd.GeoDataFrame) -> np.ndarray:
    """Compute area from polygon geometries.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame with polygon geometries.

    Returns
    -------
    np.ndarray
        Array of polygon areas.
    """
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError("Input must be a GeoDataFrame")
    return gdf.geometry.area.values


def compute_aspect_ratio_from_polygons(gdf: gpd.GeoDataFrame) -> np.ndarray:
    """Compute aspect ratio (width / height) from polygon bounding boxes.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame with polygon geometries. Must contain 'cell_id'.

    Returns
    -------
    np.ndarray
        Array of aspect ratios. NaN for multipolygons or empty geometries.
    """
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError("Input must be a GeoDataFrame")
    if "cell_id" not in gdf.columns:
        raise ValueError("'polygons' must contain a 'cell_id' column.")

    is_multi = gdf.get("is_multi", pd.Series(False, index=gdf.index))

    ar = np.full(len(gdf), np.nan)
    valid_mask = ~is_multi & ~gdf.geometry.is_empty

    for idx in gdf.index[valid_mask]:
        geom = gdf.loc[idx, "geometry"]
        bounds = geom.bounds  # (minx, miny, maxx, maxy)
        width = bounds[2] - bounds[0]
        height = bounds[3] - bounds[1]
        if height > 0:
            ar[gdf.index.get_loc(idx)] = width / height

    return ar


def compute_center_from_polygons(gdf: gpd.GeoDataFrame,
                                  coldata: pd.DataFrame) -> pd.DataFrame:
    """Compute centroid coordinates from polygon data.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame with polygon geometries.
    coldata : DataFrame
        Cell metadata to add center coordinates to.

    Returns
    -------
    DataFrame
        Updated DataFrame with center_x and center_y columns.
    """
    centroids = gdf.geometry.centroid
    cd = coldata.copy()
    cd["center_x"] = centroids.x.values
    cd["center_y"] = centroids.y.values
    return cd


def add_polygons_to_adata(adata, polygons: gpd.GeoDataFrame,
                           polygons_col: str = "polygons"):
    """Add polygon geometries to AnnData obs.

    Parameters
    ----------
    adata : AnnData
        AnnData object to add polygons to.
    polygons : GeoDataFrame
        GeoDataFrame with polygon geometries.
    polygons_col : str
        Column name for polygons in obs. Default "polygons".
    """
    import anndata as ad

    common_ids = adata.obs.index.intersection(polygons.index)
    if len(common_ids) == 0:
        raise ValueError("No matching cell IDs between adata and polygons.")
    if len(common_ids) < adata.n_obs:
        warnings.warn(
            f"{len(common_ids)}/{adata.n_obs} cells have polygon data; subsetting."
        )
    adata._inplace_subset_obs(common_ids)
    adata.obs[polygons_col] = polygons.loc[adata.obs.index, "geometry"].values
