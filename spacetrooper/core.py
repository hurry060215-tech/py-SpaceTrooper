"""SpaceTrooper class for spatial transcriptomics QC.

Provides method-chaining API wrapping the functional API.
"""

from typing import Optional, List

import anndata as ad

from .qc_metrics import spatial_per_cell_qc
from .outlier_detection import (
    compute_spatial_outlier,
    compute_threshold_flags,
    compute_outliers_qc_score,
    check_outliers,
)
from .qc_score import compute_qc_score
from .readers import read_cosmx_spe, read_merfish_spe, read_xenium_spe
from .plotting import plot_metric_hist, plot_centroids, plot_qc_score_terms


class SpaceTrooper:
    """Quality Control pipeline for image-based spatial transcriptomics.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data.
    """

    def __init__(self, adata: ad.AnnData):
        self.adata = adata

    @classmethod
    def from_cosmx(cls, dir_name: str, **kwargs) -> "SpaceTrooper":
        """Load CosMx data."""
        adata = read_cosmx_spe(dir_name, **kwargs)
        return cls(adata)

    @classmethod
    def from_merfish(cls, dir_name: str, **kwargs) -> "SpaceTrooper":
        """Load MERFISH data."""
        adata = read_merfish_spe(dir_name, **kwargs)
        return cls(adata)

    @classmethod
    def from_xenium(cls, dir_name: str, **kwargs) -> "SpaceTrooper":
        """Load Xenium data."""
        adata = read_xenium_spe(dir_name, **kwargs)
        return cls(adata)

    def spatial_per_cell_qc(self, micron_conv_fact: float = 0.12,
                            rm_zeros: bool = True,
                            neg_prob_list: Optional[List[str]] = None) -> "SpaceTrooper":
        """Compute per-cell QC metrics."""
        spatial_per_cell_qc(self.adata, micron_conv_fact=micron_conv_fact,
                           rm_zeros=rm_zeros, neg_prob_list=neg_prob_list)
        return self

    def compute_threshold_flags(self, total_threshold: float = 0,
                                ctrl_tot_ratio_threshold: float = 0.1) -> "SpaceTrooper":
        """Compute fixed-threshold flags."""
        compute_threshold_flags(self.adata, total_threshold=total_threshold,
                               ctrl_tot_ratio_threshold=ctrl_tot_ratio_threshold)
        return self

    def compute_qc_score(self, best_lambda: Optional[float] = None,
                         verbose: bool = False) -> "SpaceTrooper":
        """Compute QC score via ridge logistic regression."""
        compute_qc_score(self.adata, best_lambda=best_lambda, verbose=verbose)
        return self

    def compute_spatial_outlier(self, compute_by: str,
                                method: str = "mc",
                                mc_do_scale: bool = False,
                                scuttle_type: str = "both") -> "SpaceTrooper":
        """Compute outliers on a specific metric."""
        compute_spatial_outlier(self.adata, compute_by=compute_by,
                               method=method, mc_do_scale=mc_do_scale,
                               scuttle_type=scuttle_type)
        return self

    def plot_metric_hist(self, metric: str, **kwargs):
        """Plot histogram for a metric."""
        return plot_metric_hist(self.adata, metric=metric, **kwargs)

    def plot_centroids(self, colour_by: Optional[str] = None, **kwargs):
        """Plot spatial coordinates."""
        return plot_centroids(self.adata, colour_by=colour_by, **kwargs)

    def plot_qc_score_terms(self, **kwargs):
        """Plot QC score term maps."""
        return plot_qc_score_terms(self.adata, **kwargs)

    def __repr__(self) -> str:
        return f"SpaceTrooper({self.adata.n_obs} cells, {self.adata.n_vars} features)"
