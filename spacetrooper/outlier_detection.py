"""Outlier detection for spatial transcriptomics QC.

Ported from R/QC.R in SpaceTrooper.
Implements medcouple-based and MAD-based outlier detection.
"""

import warnings
from typing import Optional, Union, List

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import iqr as scipy_iqr


def _mc(x: np.ndarray, do_scale: bool = False) -> float:
    """Compute the medcouple — a robust measure of skewness.

    Implements the algorithm from Hubert & Vandervieren (2008).
    MC = median{ h(x_i, x_j) : x_i <= med <= x_j }
    where h(x_i, x_j) = (x_j - med + x_i - med) / (x_j - x_i)
    and for x_i == x_j: h = sign(x_i - med).

    Parameters
    ----------
    x : np.ndarray
        1D array of numeric values (NaNs removed).
    do_scale : bool
        If True, scale values before computation. Default False (as per authors).

    Returns
    -------
    float
        Medcouple value in [-1, 1].
    """
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan

    x = np.sort(x)
    med = np.median(x)

    if do_scale:
        mad = np.median(np.abs(x - med))
        if mad > 0:
            x = (x - med) / (2 * mad * 1.4826)
            med = np.median(x)

    # Split into left (<= med) and right (>= med) halves
    left = x[x <= med]
    right = x[x >= med]

    if len(left) == 0 or len(right) == 0:
        return 0.0

    # Compute h(x_i, x_j) for all pairs where x_i in left, x_j in right
    # h(x_i, x_j) = (x_i + x_j - 2*med) / (x_j - x_i)
    li = left[:, None]   # (len_left, 1)
    rj = right[None, :]  # (1, len_right)

    numerator = li + rj - 2 * med
    denominator = rj - li

    with np.errstate(divide='ignore', invalid='ignore'):
        h = np.where(denominator != 0, numerator / denominator, 0.0)

    # For x_i == x_j (when both equal med), h = sign(x_i - med) = 0
    # This is handled by the 0.0 default above

    mc_val = np.median(h)
    return float(np.clip(mc_val, -1.0, 1.0))


def _adjbox_stats(x: np.ndarray) -> dict:
    """Compute adjusted boxplot statistics (fence values).

    Uses the medcouple for skewness-adjusted fence computation.

    Parameters
    ----------
    x : np.ndarray
        1D array of numeric values.

    Returns
    -------
    dict
        Dictionary with 'fence' (lower, upper) and 'out' (outlier values).
    """
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return {"fence": (np.nan, np.nan), "out": np.array([])}

    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    mc_val = _mc(x)

    if mc_val >= 0:
        lower_fence = q1 - 1.5 * np.exp(-4 * mc_val) * iqr
        upper_fence = q3 + 1.5 * np.exp(3 * mc_val) * iqr
    else:
        lower_fence = q1 - 1.5 * np.exp(-3 * mc_val) * iqr
        upper_fence = q3 + 1.5 * np.exp(4 * mc_val) * iqr

    out = x[(x < lower_fence) | (x > upper_fence)]
    return {"fence": (lower_fence, upper_fence), "out": out}


def _is_outlier_scuttle(x: np.ndarray, type_: str = "both",
                        nmads: int = 5) -> tuple:
    """Compute MAD-based outlier detection (equivalent to scuttle::isOutlier).

    Parameters
    ----------
    x : np.ndarray
        1D array of numeric values.
    type_ : str
        One of "both", "lower", "higher".
    nmads : int
        Number of MADs from the median. Default 5.

    Returns
    -------
    tuple
        (outlier_mask, thresholds) where outlier_mask is bool array
        and thresholds is (lower, upper).
    """
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.array([], dtype=bool), (np.nan, np.nan)

    med = np.median(x)
    mad = np.median(np.abs(x - med))
    # Scale factor for normal distribution (1.4826)
    mad_scaled = mad * 1.4826

    if mad_scaled == 0:
        # Fall back to SD if MAD is 0
        mad_scaled = np.std(x, ddof=1) if len(x) > 1 else 1.0

    lower = med - nmads * mad_scaled
    upper = med + nmads * mad_scaled

    if type_ == "both":
        outlier = (x < lower) | (x > upper)
    elif type_ == "lower":
        outlier = x < lower
    elif type_ == "higher":
        outlier = x > upper
    else:
        raise ValueError(f"type_ must be 'both', 'lower', or 'higher', got '{type_}'")

    return outlier, (lower, upper)


def _check_skewness(coldata: pd.DataFrame,
                     metric_list: List[str]) -> dict:
    """Check skewness of metrics to choose outlier detection method.

    Parameters
    ----------
    coldata : DataFrame
        Cell metadata.
    metric_list : list of str
        Metrics to check.

    Returns
    -------
    dict
        Mapping of metric name to method ("mc" or "sc").
    """
    method = {}
    for metric in metric_list:
        if metric not in coldata.columns:
            continue
        vals = coldata[metric].dropna().values
        if len(vals) == 0:
            method[metric] = "sc"
            continue
        skw = stats.skew(vals)
        method[metric] = "sc" if (-1 < skw < 1) else "mc"

    # Special handling for metrics with zero-inflation
    special = {
        "log2SignalDensity": lambda cd: cd["total"] > 0,
        "log2Ctrl_total_ratio": lambda cd: cd["ctrl_total_ratio"] != 0,
    }
    for metric, mask_fn in special.items():
        if metric not in method:
            continue
        mask = mask_fn(coldata)
        vals = coldata.loc[mask, metric].dropna().values
        if len(vals) == 0:
            continue
        skw = stats.skew(vals)
        new_method = "sc" if (-1 < skw < 1) else "mc"
        if not np.isnan(skw):
            method[metric] = new_method

    return method


def compute_spatial_outlier(adata, compute_by: str,
                            method: str = "mc",
                            mc_do_scale: bool = False,
                            scuttle_type: str = "both",
                            nmads: int = 3) -> None:
    """Compute outliers based on a specified metric column.

    Adds outlier labels to adata.obs and thresholds to adata.uns.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data with QC metrics in .obs.
    compute_by : str
        Column name in adata.obs to compute outliers on.
    method : str
        One of "mc" (medcouple), "scuttle" (MAD), "both".
    mc_do_scale : bool
        Whether to scale values before medcouple computation.
    scuttle_type : str
        One of "both", "lower", "higher" for scuttle method.
    nmads : int
        Number of MADs from the median for scuttle method. Default 3.
    """
    if compute_by not in adata.obs.columns:
        raise ValueError(f"'{compute_by}' not found in adata.obs")

    cd = adata.obs
    values = cd[compute_by].values.astype(float)

    if method in ("mc", "both"):
        skw = stats.skew(values[~np.isnan(values)])
        if -1 < skw < 1:
            warnings.warn("Distribution is symmetric: mc is for asymmetric distributions. "
                          "Use scuttle instead.")

        mc_val = _mc(values, do_scale=mc_do_scale)
        if mc_val <= -0.6 or mc_val >= 0.6:
            raise ValueError(f"mc is: {mc_val:.4f}, outliers requirements not satisfied")

        stats_result = _adjbox_stats(values)
        fence = stats_result["fence"]

        outlier_labels = np.array(["NO"] * len(values), dtype=object)
        outlier_mask_low = values < fence[0]
        outlier_mask_high = values > fence[1]
        outlier_labels[outlier_mask_low] = "LOW"
        outlier_labels[outlier_mask_high] = "HIGH"

        col_name = f"{compute_by}_outlier_mc"
        adata.obs[col_name] = outlier_labels
        adata.uns[f"{col_name}_thresholds"] = np.array(fence)

    if method in ("scuttle", "both"):
        outlier_mask, thresholds = _is_outlier_scuttle(values, type_=scuttle_type, nmads=nmads)

        # Match R scuttle behavior: LOW if value <= lower fence, HIGH otherwise
        outlier_labels = np.array(["NO"] * len(values), dtype=object)
        outlier_labels[outlier_mask & (values <= thresholds[0])] = "LOW"
        outlier_labels[outlier_mask & ~(values <= thresholds[0])] = "HIGH"

        col_name = f"{compute_by}_outlier_sc"
        adata.obs[col_name] = outlier_labels
        adata.uns[f"{col_name}_thresholds"] = np.array(thresholds)


def compute_threshold_flags(adata, total_threshold: float = 0,
                            ctrl_tot_ratio_threshold: float = 0.1) -> None:
    """Compute flagged cells using fixed thresholds.

    Adds boolean flags to adata.obs.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data.
    total_threshold : float
        Threshold for total counts. Default 0.
    ctrl_tot_ratio_threshold : float
        Threshold for control-to-total ratio. Default 0.1.
    """
    if "total" not in adata.obs.columns:
        raise ValueError("'total' not found in adata.obs")
    if "ctrl_total_ratio" not in adata.obs.columns:
        raise ValueError("'ctrl_total_ratio' not found in adata.obs")

    adata.obs["is_zero_counts"] = adata.obs["total"] == total_threshold
    adata.obs["is_ctrl_tot_outlier"] = (
        adata.obs["ctrl_total_ratio"] > ctrl_tot_ratio_threshold
    )
    adata.obs["threshold_flags"] = (
        adata.obs["is_ctrl_tot_outlier"] & adata.obs["is_zero_counts"]
    )


def compute_outliers_qc_score(adata, metric_list: Optional[List[str]] = None) -> None:
    """Compute outlier cells for each QC score formula variable.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data with QC metrics in .obs.
    metric_list : list of str, optional
        Metrics to compute outliers for. Default: log2SignalDensity,
        Area_um, log2AspectRatio, log2Ctrl_total_ratio.
    """
    if metric_list is None:
        metric_list = ["log2SignalDensity", "Area_um",
                        "log2AspectRatio", "log2Ctrl_total_ratio"]

    cd = adata.obs
    raw_method = _check_skewness(cd, metric_list)
    # Map shorthand: "sc" -> "scuttle", "mc" -> "mc"
    method = {k: ("scuttle" if v == "sc" else v) for k, v in raw_method.items()}
    # Also keep the short form for column name lookups
    method_short = {k: v for k, v in raw_method.items()}

    # Special handling for zero-inflated metrics
    spec = {
        "log2SignalDensity": {
            "idx": adata.obs["total"] > 0,
            "zero": adata.obs["total"] == 0,
            "tweak_lower": True,
        },
        "log2Ctrl_total_ratio": {
            "idx": adata.obs["ctrl_total_ratio"] != 0,
            "zero": adata.obs["ctrl_total_ratio"] == 0,
            "tweak_lower": False,
        },
    }

    # Remove metrics with None/NA method
    na_metrics = [m for m in method if method.get(m) is None]
    for m in na_metrics:
        warnings.warn(f"Metric {m} produced NA skewness values. Removing from formula.")
        del method[m]

    # Process zero-inflated metrics
    for var, s in spec.items():
        if var not in method:
            continue
        adata_sub = adata[s["idx"]].copy()
        compute_spatial_outlier(adata_sub, compute_by=var, method=method[var])
        out_var = f"{var}_outlier_{method_short[var]}"

        # Get fences
        thresholds_key = f"{out_var}_thresholds"
        thresholds = adata_sub.uns[thresholds_key]
        low_thr = thresholds[0]
        high_thr = thresholds[1]

        if s["tweak_lower"]:
            if low_thr < adata_sub.obs[var].min():
                low_thr = np.percentile(adata.obs[var].dropna().values, 1)

        # Apply to all cells
        train_col = f"{var}_outlier_train"
        labels = np.where(s["zero"], "NO",
                         np.where(adata.obs[var] < low_thr, "LOW",
                                 np.where(adata.obs[var] > high_thr, "HIGH", "NO")))
        adata.obs[train_col] = labels

        # Store thresholds
        thr = thresholds.copy()
        if s["tweak_lower"]:
            thr[0] = low_thr
        adata.uns[f"{train_col}_thresholds"] = thr

    # Process non-zero-inflated metrics
    sub_method = {k: v for k, v in method.items()
                  if k not in ("log2SignalDensity", "log2Ctrl_total_ratio")}
    for var in sub_method:
        compute_spatial_outlier(adata, compute_by=var, method=sub_method[var])

    # Build formula_variables mapping
    out_var = {}
    for m in method:
        if m == "log2SignalDensity":
            out_var[m] = "log2SignalDensity_outlier_train"
        elif m == "log2Ctrl_total_ratio":
            out_var[m] = "log2Ctrl_total_ratio_outlier_train"
        else:
            out_var[m] = f"{m}_outlier_{method_short[m]}"

    adata.uns["formula_variables"] = out_var


def check_outliers(adata, verbose: bool = False) -> None:
    """Check if computed outliers meet minimum numerical requirements.

    Removes variables from formula if outlier count < 0.1% of cells.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data.
    verbose : bool
        Print outlier counts.
    """
    out_var = dict(adata.uns.get("formula_variables", {}))
    cd = adata.obs

    if verbose:
        for var in out_var:
            col = out_var[var]
            if col in cd.columns:
                print(f"Outliers found for {var}:")
                print(cd[col].value_counts().to_string())

    if "log2SignalDensity" not in out_var:
        raise ValueError("log2SignalDensity is not in the QC score formula. "
                         "QC score cannot be computed.")

    n_min = len(adata) * 0.001

    config = {
        "log2SignalDensity": {
            "pattern": "log2SignalDensity_outlier",
            "remove": "log2SignalDensity_outlier_train",
            "label": "LOW",
            "action": "stop",
        },
        "Area_um": {
            "pattern": "Area_um_outlier",
            "remove": "Area_um_outlier",
            "label": "HIGH",
            "action": "warn",
        },
        "log2Ctrl_total_ratio": {
            "pattern": "log2Ctrl_total_ratio_outlier",
            "remove": "log2Ctrl_total_ratio_outlier_train",
            "label": "HIGH",
            "action": "warn",
        },
    }

    for var, cfg in config.items():
        if var not in out_var:
            continue
        col = out_var[var]
        if col not in cd.columns:
            continue

        counts = cd[col].value_counts()
        count = counts.get(cfg["label"], 0)

        if count < n_min:
            msg = f"Not enough outlier cells for {var}.\n"
            if cfg["action"] == "stop":
                msg += "QC score cannot be computed without log2SignalDensity outliers."
                raise ValueError(msg)
            else:
                msg += "This variable will not be used in the final formula."
                warnings.warn(msg)
                if cfg["remove"] in out_var:
                    del out_var[cfg["remove"]]

    # Handle log2AspectRatio for CosMx
    var = "log2AspectRatio"
    is_cosmx = adata.uns.get("technology", "") in (
        "Nanostring_CosMx", "Nanostring_CosMx_Protein"
    )
    if is_cosmx and var in out_var:
        col = out_var[var]
        if col in cd.columns:
            counts = cd[col].value_counts()
            low = counts.get("LOW", 0)
            high = counts.get("HIGH", 0)
            if low < n_min and high < n_min:
                warnings.warn(f"Not enough outlier cells for {var}.\n"
                              "This variable will not be used in the final formula.")
                out_var = {k: v for k, v in out_var.items()
                          if not v.startswith("log2AspectRatio_outlier")}
    else:
        out_var = {k: v for k, v in out_var.items()
                  if not v.startswith("log2AspectRatio_outlier")}

    adata.uns["formula_variables"] = out_var
