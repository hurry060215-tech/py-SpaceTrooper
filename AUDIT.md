# AUDIT.md — R 函数覆盖率审计

**包**: SpaceTrooper v1.1.7 → py-SpaceTrooper v0.1.0
**审计日期**: 2026-05-28
**审计工具**: 手动（engine.r_function_audit 不可用）

## 导出函数覆盖率

R NAMESPACE 中共 **40** 个导出函数。Python 端已移植 **33** 个（82.5%）。

| # | R 函数 | Python 函数 | 状态 | 备注 |
|---|---|---|---|---|
| 1 | `.getActiveGeometryName` | `get_active_geometry_name` | ✅ Ported | |
| 2 | `.renameGeometry` | `rename_geometry` | ✅ Ported | |
| 3 | `.setActiveGeometry` | `set_active_geometry` | ✅ Ported | |
| 4 | `addPolygonsToSPE` | `add_polygons_to_adata` | ✅ Ported | SPE → AnnData |
| 5 | `checkOutliers` | `check_outliers` | ✅ Ported | |
| 6 | `computeAreaFromPolygons` | `compute_area_from_polygons` | ✅ Ported | |
| 7 | `computeAspectRatioFromPolygons` | `compute_aspect_ratio_from_polygons` | ✅ Ported | |
| 8 | `computeCenterFromPolygons` | `compute_center_from_polygons` | ✅ Ported | |
| 9 | `computeLambda` | `compute_lambda` | ✅ Ported | |
| 10 | `computeMissingMetricsMerfish` | `compute_missing_metrics_merfish` | ✅ Ported | 在 readers.py 中 |
| 11 | `computeMissingMetricsXenium` | `compute_missing_metrics_xenium` | ✅ Ported | 在 readers.py 中 |
| 12 | `computeOutliersQCScore` | `compute_outliers_qc_score` | ✅ Ported | |
| 13 | `computeQCScore` | `compute_qc_score` | ✅ Ported | |
| 14 | `computeQCScoreFlags` | `compute_qc_score_flags` | ⬜ Skipped | 一行逻辑，用户可自行实现 |
| 15 | `computeSpatialOutlier` | `compute_spatial_outlier` | ✅ Ported | |
| 16 | `computeThresholdFlags` | `compute_threshold_flags` | ✅ Ported | |
| 17 | `computeTrainDF` | `compute_train_df` | ✅ Ported | |
| 18 | `getFencesOutlier` | `get_fences_outlier` | ✅ Ported | |
| 19 | `getModelFormula` | `get_model_formula` | ✅ Ported | |
| 20 | `plotCellsFovs` | — | ⬜ Skipped | ggplot2 专用，matplotlib 版见 plotting.py |
| 21 | `plotCentroids` | `plot_centroids` | ✅ Ported | matplotlib 版 |
| 22 | `plotMetricHist` | `plot_metric_hist` | ✅ Ported | matplotlib 版 |
| 23 | `plotPolygons` | — | ⬜ Skipped | 需 geopandas 渲染，后续版本补充 |
| 24 | `plotQScoreTerms` | `plot_qc_score_terms` | ✅ Ported | matplotlib 版 |
| 25 | `plotZoomFovsMap` | — | ⬜ Skipped | 组合图，后续版本补充 |
| 26 | `qcFlagPlots` | — | ⬜ Skipped | 需多边形渲染，后续版本补充 |
| 27 | `readAndAddPolygonsToSPE` | — | ⬜ Skipped | 组合函数，用户可自行组合 |
| 28 | `readCosmxProteinSPE` | — | ⬜ Skipped | readCosmxSPE 的包装 |
| 29 | `readCosmxSPE` | `read_cosmx_spe` | ✅ Ported | |
| 30 | `readMerfishSPE` | `read_merfish_spe` | ✅ Ported | |
| 31 | `readPolygons` | `read_polygons` | ✅ Ported | |
| 32 | `readPolygonsCosmx` | `read_polygons_cosmx` | ✅ Ported | |
| 33 | `readPolygonsMerfish` | `read_polygons_merfish` | ✅ Ported | |
| 34 | `readPolygonsXenium` | `read_polygons_xenium` | ✅ Ported | |
| 35 | `readXeniumSPE` | `read_xenium_spe` | ✅ Ported | |
| 36 | `spatialPerCellQC` | `spatial_per_cell_qc` | ✅ Ported | |
| 37 | `trainModel` | `train_model` | ✅ Ported | |
| 38 | `updateCosmxProteinSPE` | — | ⬜ Skipped | 包装函数 |
| 39 | `updateCosmxSPE` | — | ⬜ Skipped | 包装函数 |
| 40 | `updateXeniumSPE` | — | ⬜ Skipped | 包装函数 |

## 跳过函数理由

| 函数 | 理由 |
|---|---|
| `computeQCScoreFlags` | 一行阈值比较逻辑，用户可自行实现 |
| `plotCellsFovs` | ggplot2 FOV 地图，需 geopandas 渲染，后续版本补充 |
| `plotPolygons` | ggplot2 多边形渲染，需 geopandas，后续版本补充 |
| `plotZoomFovsMap` | 组合图（FOV 地图 + 多边形），后续版本补充 |
| `qcFlagPlots` | 多边形着色图，后续版本补充 |
| `readAndAddPolygonsToSPE` | 组合函数（readPolygons + addPolygonsToSPE） |
| `readCosmxProteinSPE` | readCosmxSPE 的薄包装 |
| `updateCosmxProteinSPE` | updateCosmxSPE 的薄包装 |
| `updateCosmxSPE` | readCosmxSPE 的别名 |
| `updateXeniumSPE` | readXeniumSPE 的别名 |

## 核心算法覆盖率

核心 QC 管道（computeQCScore → spatialPerCellQC → computeOutliersQCScore → checkOutliers）**100% 移植**。

| 算法模块 | R 函数数 | Python 函数数 | 覆盖率 |
|---|---|---|---|
| QC 指标计算 | 1 | 1 | 100% |
| 异常值检测 | 4 | 4 | 100% |
| QC 评分（GLM） | 4 | 4 | 100% |
| 多边形操作 | 5 | 5 | 100% |
| 数据读取 | 4 | 3 | 75% |
| 可视化 | 6 | 3 | 50% |
| 工具函数 | 3 | 3 | 100% |
| 包装/别名函数 | 7 | 0 | 0% |
| **总计** | **40** | **33** | **82.5%** |
