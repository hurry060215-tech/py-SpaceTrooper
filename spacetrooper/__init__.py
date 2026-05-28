"""SpaceTrooper — Quality Control for image-based spatial transcriptomics.

Python port of the R Bioconductor package SpaceTrooper.
"""

from .qc_metrics import spatial_per_cell_qc
from .outlier_detection import (
    compute_spatial_outlier,
    compute_threshold_flags,
    compute_outliers_qc_score,
    check_outliers,
)
from .utils import get_fences_outlier
from .qc_score import (
    compute_qc_score,
    compute_train_df,
    compute_lambda,
    train_model,
    get_model_formula,
)
from .utils import (
    get_active_geometry_name,
    set_active_geometry,
    rename_geometry,
)
from .polygons import (
    read_polygons,
    read_polygons_cosmx,
    read_polygons_xenium,
    read_polygons_merfish,
    add_polygons_to_adata,
    compute_area_from_polygons,
    compute_aspect_ratio_from_polygons,
    compute_center_from_polygons,
)
from .readers import read_cosmx_spe, read_merfish_spe, read_xenium_spe
from .core import SpaceTrooper
from .plotting import plot_metric_hist, plot_centroids, plot_qc_score_terms

__all__ = [
    # Class API
    "SpaceTrooper",
    # QC metrics
    "spatial_per_cell_qc",
    # Outlier detection
    "compute_spatial_outlier",
    "compute_threshold_flags",
    "compute_outliers_qc_score",
    "check_outliers",
    "get_fences_outlier",
    # QC score
    "compute_qc_score",
    "compute_train_df",
    "compute_lambda",
    "train_model",
    "get_model_formula",
    # Utils
    "get_active_geometry_name",
    "set_active_geometry",
    "rename_geometry",
    # Polygons
    "read_polygons",
    "read_polygons_cosmx",
    "read_polygons_xenium",
    "read_polygons_merfish",
    "add_polygons_to_adata",
    "compute_area_from_polygons",
    "compute_aspect_ratio_from_polygons",
    "compute_center_from_polygons",
    # Readers
    "read_cosmx_spe",
    "read_merfish_spe",
    "read_xenium_spe",
    # Plotting
    "plot_metric_hist",
    "plot_centroids",
    "plot_qc_score_terms",
]
__version__ = "0.1.1"
