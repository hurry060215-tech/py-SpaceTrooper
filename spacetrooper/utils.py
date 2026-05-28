"""Utility functions for geometry management and outlier fence extraction.

Ported from R/geometry.R and R/utils.R in SpaceTrooper.
"""

import numpy as np
import pandas as pd
import geopandas as gpd


def get_active_geometry_name(gdf: gpd.GeoDataFrame) -> str:
    """Return the name of the active geometry column.

    Parameters
    ----------
    gdf : GeoDataFrame
        A GeoDataFrame with an active geometry column.

    Returns
    -------
    str
        Name of the active geometry column.
    """
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError("Input must be a GeoDataFrame")
    return gdf.geometry.name


def set_active_geometry(gdf: gpd.GeoDataFrame, name: str) -> gpd.GeoDataFrame:
    """Set the active geometry column of a GeoDataFrame.

    Parameters
    ----------
    gdf : GeoDataFrame
        A GeoDataFrame with multiple geometry columns.
    name : str
        Name of the geometry column to activate.

    Returns
    -------
    GeoDataFrame
        GeoDataFrame with the specified geometry set as active.
    """
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError("Input must be a GeoDataFrame")
    if name not in gdf.columns:
        raise ValueError(f"Geometry '{name}' not found in GeoDataFrame columns")
    gdf = gdf.set_geometry(name)
    return gdf


def rename_geometry(gdf: gpd.GeoDataFrame, from_name: str, to_name: str,
                    activate: bool = False) -> gpd.GeoDataFrame:
    """Rename a geometry column in a GeoDataFrame.

    Parameters
    ----------
    gdf : GeoDataFrame
        A GeoDataFrame containing the geometry to rename.
    from_name : str
        Current name of the geometry column.
    to_name : str
        New name for the geometry column.
    activate : bool
        If True, set the renamed geometry as active. Default False.

    Returns
    -------
    GeoDataFrame
        GeoDataFrame with the renamed geometry.
    """
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError("Input must be a GeoDataFrame")
    if from_name not in gdf.columns:
        raise ValueError(f"Geometry '{from_name}' not found in GeoDataFrame columns")

    was_active = (gdf.geometry.name == from_name)
    gdf = gdf.rename(columns={from_name: to_name})

    if was_active or activate:
        gdf = gdf.set_geometry(to_name)

    return gdf


def get_fences_outlier(adata, fences_of: str, high_low: str = "both",
                       decimal_round: int = None) -> np.ndarray:
    """Extract threshold (fence) values from outlier detection results.

    Parameters
    ----------
    adata : AnnData
        An AnnData object with outlier results in .uns.
    fences_of : str
        Name of the column in adata.obs containing outlier labels.
    high_low : str
        One of "both", "lower", "higher". Default "both".
    decimal_round : int, optional
        Number of decimal places for rounding.

    Returns
    -------
    np.ndarray or float
        Fence values (lower and upper thresholds).
    """
    thresholds_key = f"{fences_of}_thresholds"
    if thresholds_key not in adata.uns:
        raise KeyError(f"Thresholds '{thresholds_key}' not found in adata.uns")

    fences = np.array(adata.uns[thresholds_key], dtype=float)
    if decimal_round is not None:
        fences = np.round(fences, decimal_round)

    if high_low == "both":
        return fences
    elif high_low == "lower":
        return fences[0]
    elif high_low == "higher":
        return fences[1]
    else:
        raise ValueError(f"high_low must be 'both', 'lower', or 'higher', got '{high_low}'")
