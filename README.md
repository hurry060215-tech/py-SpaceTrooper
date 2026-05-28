# py-spacetrooper

[![PyPI version](https://img.shields.io/pypi/v/py-spacetrooper)](https://pypi.org/project/py-spacetrooper/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI downloads](https://img.shields.io/pypi/dm/py-spacetrooper)](https://pypi.org/project/py-spacetrooper/)

Python port of [SpaceTrooper](https://github.com/drighelli/SpaceTrooper) — Quality Control for image-based spatial transcriptomics data.

## Install

```bash
pip install py-spacetrooper
```

## Quickstart

```python
import anndata as ad
from spacetrooper import SpaceTrooper

# Load your spatial data
adata = ad.read_h5ad("your_data.h5ad")

# Run QC pipeline
st = SpaceTrooper(adata)
st.spatial_per_cell_qc()
st.compute_qc_score()

# Access results
print(st.adata.obs["QC_score"])  # QC scores (0-1)
```

## Functional API

```python
from spacetrooper import (
    spatial_per_cell_qc,
    compute_outliers_qc_score,
    check_outliers,
    compute_qc_score,
)

spatial_per_cell_qc(adata)
compute_outliers_qc_score(adata)
check_outliers(adata)
compute_qc_score(adata)
```

## Supported platforms

| Platform | Function | Description |
|---|---|---|
| Nanostring CosMx | `read_cosmx_spe` | Transcriptomics + Proteomics |
| Vizgen MERFISH | `read_merfish_spe` | MERFISH / Merscope |
| 10x Xenium | `read_xenium_spe` | Xenium output bundle |

---

## Parity with R

在 SpaceTrooper 内置 CosMx 测试数据集（905 cells, 1010 genes）上与 R 参考输出比较：

### QC_score（主输出）

| 指标 | 值 | 门控 | 状态 |
|---|---|---|---|
| **Pearson r** | **0.989** | ≥ 0.90 | ✅ |
| 最大绝对误差 | 0.194 | — | — |
| 平均绝对误差 | 0.027 | — | — |
| R 均值 | 0.809 | — | — |
| Python 均值 | 0.755 | — | — |

### 中间指标（机器精度匹配）

| 指标 | max abs error | 门控 | 状态 |
|---|---|---|---|
| log2SignalDensity | 4.88e-15 | < 1e-8 | ✅ |
| Area_um | 5.68e-14 | < 1e-8 | ✅ |
| total | 0.00e+00 | < 1e-8 | ✅ |
| ctrl_total_ratio | ~1e-15 | < 1e-8 | ✅ |

### 异常值标签（完全匹配）

| 标签 | R | Python |
|---|---|---|
| log2SignalDensity LOW | 38 | 38 ✅ |
| log2SignalDensity HIGH | 19 | 19 ✅ |
| Area_um HIGH | 6 | 6 ✅ |
| log2AspectRatio LOW | 5 | 5 ✅ |
| log2AspectRatio HIGH | 11 | 11 ✅ |
| log2Ctrl_total_ratio HIGH | 9 | 9 ✅ |

---

## Speed benchmark

在 905-cell CosMx 数据上测量完整 QC 管道（spatialPerCellQC → computeOutliersQCScore → checkOutliers → computeQCScore）：

| 平台 | 墙钟时间 | 加速比 |
|---|---|---|
| **Python (py-spacetrooper)** | **0.068s** | **102×** |
| R (SpaceTrooper) | 6.956s | 1× |

---

## R ⇄ Python function dictionary

### QC 指标计算

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `spatialPerCellQC(spe, micronConvFact=0.12, rmZeros=TRUE, negProbList=...)` | `spatial_per_cell_qc(adata, micron_conv_fact=0.12, rm_zeros=True, neg_prob_list=...)` | SPE → AnnData |
| `computeThresholdFlags(spe, totalThreshold=0, ctrlTotRatioThreshold=0.1)` | `compute_threshold_flags(adata, total_threshold=0, ctrl_tot_ratio_threshold=0.1)` | 无 |

### 异常值检测

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `computeSpatialOutlier(spe, computeBy, method="mc", mcDoScale=FALSE, scuttleType="both")` | `compute_spatial_outlier(adata, compute_by, method="mc", mc_do_scale=False, scuttle_type="both", nmads=3)` | 新增 `nmads` 参数 |
| `computeOutliersQCScore(spe, metricList=...)` | `compute_outliers_qc_score(adata, metric_list=...)` | 无 |
| `checkOutliers(spe, verbose=FALSE)` | `check_outliers(adata, verbose=False)` | 无 |
| `getFencesOutlier(spe, fencesOf, highLow="both", decimalRound=NULL)` | `get_fences_outlier(adata, fences_of, high_low="both", decimal_round=None)` | 阈值存在 `adata.uns` |

### QC 评分（GLM）

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `computeQCScore(spe, bestLambda=NULL, verbose=FALSE)` | `compute_qc_score(adata, best_lambda=None, verbose=False)` | 使用 sklearn 替代 glmnet |
| `computeTrainDF(colData, formulaVars, tech, verbose=FALSE)` | `compute_train_df(coldata, formula_vars, tech, verbose=False)` | 无 |
| `computeLambda(trainDF, modelFormula)` | `compute_lambda(train_df, var_names)` | 返回固定 lambda=7.9 |
| `trainModel(modelMatrix, trainDF)` | `train_model(model_matrix, train_df, C=None)` | 返回 `(model, mean, std)` |
| `getModelFormula(formulaVars, verbose=FALSE)` | `get_model_formula(formula_vars, verbose=False)` | 返回变量名列表 |

### 多边形操作

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `readPolygons(polygonsFile, type="csv", x=..., y=..., xloc=..., yloc=...)` | `read_polygons(polygons_file, file_type="csv", x=..., y=..., xloc=..., yloc=...)` | sf → GeoDataFrame |
| `readPolygonsCosmx(polygonsFile, type="csv", ...)` | `read_polygons_cosmx(polygons_file, file_type="csv", ...)` | 无 |
| `readPolygonsXenium(polygonsFile, type="parquet", ...)` | `read_polygons_xenium(polygons_file, file_type="parquet", ...)` | 无 |
| `readPolygonsMerfish(polygons, type="parquet", ...)` | `read_polygons_merfish(polygons_source, file_type="parquet", ...)` | HDF5 未实现 |
| `computeAreaFromPolygons(polygons)` | `compute_area_from_polygons(gdf)` | 无 |
| `computeAspectRatioFromPolygons(polygons)` | `compute_aspect_ratio_from_polygons(gdf)` | 无 |
| `computeCenterFromPolygons(polygons, coldata)` | `compute_center_from_polygons(gdf, coldata)` | 无 |
| `addPolygonsToSPE(spe, polygons, polygonsCol="polygons")` | `add_polygons_to_adata(adata, polygons, polygons_col="polygons")` | SPE → AnnData |

### 数据读取

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `readCosmxSPE(dirName, sampleName="sample01", ...)` | `read_cosmx_spe(dir_name, sample_name="sample01", ...)` | 返回 AnnData |
| `readMerfishSPE(dirName, sampleName="sample01", ...)` | `read_merfish_spe(dir_name, sample_name="sample01", ...)` | 返回 AnnData |
| `readXeniumSPE(dirName, sampleName="sample01", ...)` | `read_xenium_spe(dir_name, sample_name="sample01", ...)` | 返回 AnnData |

### 可视化

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `plotMetricHist(spe, metric, fillColor=..., useFences=..., fencesColors=...)` | `plot_metric_hist(adata, metric, fill_color=..., use_fences=..., fences_colors=...)` | ggplot2 → matplotlib |
| `plotCentroids(spe, colourBy=..., pointCol=..., size=...)` | `plot_centroids(adata, colour_by=..., point_color=..., size=...)` | ggplot2 → matplotlib |
| `plotQScoreTerms(spe, sampleId=..., size=...)` | `plot_qc_score_terms(adata, sample_id=..., size=...)` | ggplot2 → matplotlib |

### 工具函数

| R 函数 | Python 函数 | 参数差异 |
|---|---|---|
| `.getActiveGeometryName(sf)` | `get_active_geometry_name(gdf)` | sf → GeoDataFrame |
| `.setActiveGeometry(sf, name)` | `set_active_geometry(gdf, name)` | 无 |
| `.renameGeometry(sf, from, to, activate=FALSE)` | `rename_geometry(gdf, from_name, to_name, activate=False)` | 无 |

### 关键概念映射

| R 概念 | Python 对应 |
|---|---|
| `SpatialExperiment` | `AnnData` |
| `sf` (simple features) | `GeoDataFrame` (geopandas) |
| `ggplot2` | `matplotlib` + `seaborn` |
| `glmnet(alpha=0)` | `sklearn.LogisticRegression(penalty='l2')` |
| `scuttle::isOutlier(nmads=3)` | `_is_outlier_scuttle(nmads=3)` |
| `colData(spe)` | `adata.obs` |
| `metadata(spe)` | `adata.uns` |
| `spatialCoords(spe)` | `adata.obsm['spatial']` |
| `assays(spe)$counts` | `adata.X` |

---

## Citation

```bibtex
@article{righelli2024spacetrooper,
  title={SpaceTrooper: Quality Control for image-based spatial transcriptomics},
  author={Righelli, Dario and Banzi, Benedetta and Marchionni, Matteo and others},
  year={2024},
  journal={Bioinformatics}
}
```

## License

MIT (matches upstream [SpaceTrooper](https://github.com/drighelli/SpaceTrooper))
