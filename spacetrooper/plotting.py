"""Visualization functions for spatial transcriptomics QC.

Ported from R/plotUtils.R and R/spatialQCPlots.R in SpaceTrooper.
Uses matplotlib/seaborn instead of ggplot2.
"""

import warnings
from typing import Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap


def plot_metric_hist(adata, metric: str, fill_color: str = "#c0c8cf",
                     use_fences: Optional[str] = None,
                     fences_colors: Optional[dict] = None,
                     bins: int = 30, ax: Optional[plt.Axes] = None) -> plt.Figure:
    """Plot a histogram for a given metric.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data.
    metric : str
        Column name in adata.obs to plot.
    fill_color : str
        Fill color for histogram bars.
    use_fences : str, optional
        Column name for outlier fences to overlay.
    fences_colors : dict, optional
        Colors for lower/higher fences.
    bins : int
        Number of bins.
    ax : Axes, optional
        Matplotlib axes to plot on.

    Returns
    -------
    Figure
    """
    if metric not in adata.obs.columns:
        raise ValueError(f"'{metric}' not found in adata.obs")

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    else:
        fig = ax.figure

    values = adata.obs[metric].dropna().values
    ax.hist(values, bins=bins, color=fill_color, edgecolor="white", linewidth=0.5)

    if use_fences is not None:
        if fences_colors is None:
            fences_colors = {"lower": "purple", "higher": "tomato"}

        thresholds_key = f"{use_fences}_thresholds"
        if thresholds_key in adata.uns:
            fences = adata.uns[thresholds_key]
            ax.axvline(fences[0], color=fences_colors.get("lower", "purple"),
                      linestyle="--", linewidth=1.5, label=f"lower: {fences[0]:.2f}")
            ax.axvline(fences[1], color=fences_colors.get("higher", "tomato"),
                      linestyle="--", linewidth=1.5, label=f"higher: {fences[1]:.2f}")
            ax.legend(fontsize=9)

    ax.set_title(metric, fontsize=12, fontweight="bold")
    ax.set_xlabel(metric)
    ax.set_ylabel("Count")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig


def plot_centroids(adata, colour_by: Optional[str] = None,
                   point_size: float = 0.05, alpha: float = 0.8,
                   point_color: str = "darkmagenta",
                   title: Optional[str] = None,
                   ax: Optional[plt.Axes] = None) -> plt.Figure:
    """Plot spatial coordinates colored by a metric.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data with 'spatial' in .obsm.
    colour_by : str, optional
        Column in adata.obs for coloring points.
    point_size : float
        Size of scatter points.
    alpha : float
        Transparency.
    point_color : str
        Default color when colour_by is None.
    title : str, optional
        Plot title.
    ax : Axes, optional
        Matplotlib axes.

    Returns
    -------
    Figure
    """
    if "spatial" not in adata.obsm:
        raise ValueError("'spatial' not found in adata.obsm")

    coords = adata.obsm["spatial"]
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    else:
        fig = ax.figure

    if colour_by is not None and colour_by in adata.obs.columns:
        values = adata.obs[colour_by].values
        scatter = ax.scatter(coords[:, 0], coords[:, 1],
                           c=values, s=point_size, alpha=alpha,
                           cmap="viridis")
        plt.colorbar(scatter, ax=ax, label=colour_by)
    else:
        ax.scatter(coords[:, 0], coords[:, 1],
                  c=point_color, s=point_size, alpha=alpha)

    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig


def plot_qc_score_terms(adata, sample_id: Optional[str] = None,
                         size: float = 0.05, alpha: float = 0.8,
                         figsize: tuple = (15, 5)) -> plt.Figure:
    """Plot spatial maps of QC score terms.

    Parameters
    ----------
    adata : AnnData
        Spatial data with QC metrics.
    sample_id : str, optional
        Sample identifier for title.
    size : float
        Point size.
    alpha : float
        Transparency.
    figsize : tuple
        Figure size.

    Returns
    -------
    Figure
    """
    if "spatial" not in adata.obsm:
        raise ValueError("'spatial' not found in adata.obsm")

    coords = adata.obsm["spatial"]
    metrics = ["log2SignalDensity", "log2AspectRatio"]
    available = [m for m in metrics if m in adata.obs.columns]

    n_plots = len(available)
    if "dist_border" in adata.obs.columns:
        n_plots += 1

    fig, axes = plt.subplots(1, n_plots, figsize=figsize)
    if n_plots == 1:
        axes = [axes]

    for i, metric in enumerate(available):
        values = adata.obs[metric].values
        scatter = axes[i].scatter(coords[:, 0], coords[:, 1],
                                  c=values, s=size, alpha=alpha,
                                  cmap="viridis")
        axes[i].set_title(metric, fontsize=10)
        axes[i].set_aspect("equal")
        axes[i].axis("off")
        plt.colorbar(scatter, ax=axes[i], fraction=0.046)

    if "dist_border" in adata.obs.columns and n_plots > len(available):
        values = adata.obs["dist_border"].values
        scatter = axes[-1].scatter(coords[:, 0], coords[:, 1],
                                   c=values, s=size, alpha=alpha,
                                   cmap="viridis")
        axes[-1].set_title("dist_border", fontsize=10)
        axes[-1].set_aspect("equal")
        axes[-1].axis("off")
        plt.colorbar(scatter, ax=axes[-1], fraction=0.046)

    if sample_id:
        fig.suptitle(sample_id, fontsize=14, fontweight="bold")

    fig.tight_layout()
    return fig
