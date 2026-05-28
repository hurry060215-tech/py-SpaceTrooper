#!/usr/bin/env Rscript
# Simplified R reference driver for SpaceTrooper parity testing.
# Constructs SpatialExperiment manually from CSV files.
#
# Usage: Rscript r_reference_driver.R <cosmx_dir> <output_json>

.libPaths(unique(c("C:/Users/17904/Documents/R/win-library/4.5", .libPaths())))

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(robustbase)
    library(e1071)
    library(scuttle)
    library(glmnet)
    library(dplyr)
    library(jsonlite)
})

# Source SpaceTrooper functions directly
st_path <- "D:/test/SpaceTrooper-devel/SpaceTrooper-devel/R"
for (f in list.files(st_path, pattern="\\.R$", full.names=TRUE)) {
    tryCatch(source(f), error=function(e) NULL)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript r_reference_driver.R <cosmx_dir> <output_json>")
}

cosmx_dir <- args[1]
output_path <- args[2]

cat("Reading CosMx data from:", cosmx_dir, "\n")

# Read CSV files
expr_file <- list.files(cosmx_dir, pattern="exprMat_file\\.csv", full.names=TRUE)
meta_file <- list.files(cosmx_dir, pattern="metadata_file\\.csv", full.names=TRUE)
fov_file  <- list.files(cosmx_dir, pattern="fov_positions_file\\.csv", full.names=TRUE)

counts_df <- read.csv(expr_file, check.names=FALSE)
meta_df   <- read.csv(meta_file, check.names=FALSE)
fov_pos   <- read.csv(fov_file, check.names=FALSE)

# Build cell IDs
counts_df$cell_id <- paste0("f", counts_df$fov, "_c", counts_df$cell_ID)
meta_df$cell_id   <- paste0("f", meta_df$fov, "_c", meta_df$cell_ID)

# Align to common cells
common_ids <- intersect(counts_df$cell_id, meta_df$cell_id)
counts_df <- counts_df[match(common_ids, counts_df$cell_id), ]
meta_df   <- meta_df[match(common_ids, meta_df$cell_id), ]

# Extract count matrix
feature_cols <- setdiff(colnames(counts_df), c("fov", "cell_ID", "cell_id"))
count_matrix <- t(as.matrix(counts_df[, feature_cols]))
colnames(count_matrix) <- common_ids

# Build colData
rownames(meta_df) <- meta_df$cell_id

# Create SpatialExperiment
spe <- SpatialExperiment(
    assay=list(counts=count_matrix),
    colData=DataFrame(meta_df),
    spatialCoordsNames=c("CenterX_global_px", "CenterY_global_px"),
    sample_id="DBKero_Tiny"
)

# Set metadata
metadata(spe) <- list(
    fov_positions=fov_pos,
    fov_dim=c(xdim=4256, ydim=4256),
    polygons=list.files(cosmx_dir, pattern="polygons\\.csv", full.names=TRUE),
    technology="Nanostring_CosMx"
)

colnames(spe) <- common_ids
spe$cell_id <- common_ids

cat("SPE created:", ncol(spe), "cells,", nrow(spe), "genes\n")

# Run QC pipeline
set.seed(42)

cat("Running spatialPerCellQC...\n")
spe <- spatialPerCellQC(spe)

cat("Running computeOutliersQCScore...\n")
spe <- computeOutliersQCScore(spe)

cat("Running checkOutliers...\n")
spe <- checkOutliers(spe, verbose=TRUE)

cat("Running computeQCScore...\n")
spe <- computeQCScore(spe, verbose=TRUE)

# Extract results
result <- list(
    QC_score = as.numeric(colData(spe)$QC_score),
    cell_id = as.character(colData(spe)$cell_id),
    log2SignalDensity = as.numeric(colData(spe)$log2SignalDensity),
    Area_um = as.numeric(colData(spe)$Area_um),
    ctrl_total_ratio = as.numeric(colData(spe)$ctrl_total_ratio),
    log2Ctrl_total_ratio = as.numeric(colData(spe)$log2Ctrl_total_ratio),
    log2AspectRatio = as.numeric(colData(spe)$log2AspectRatio),
    signal_density = as.numeric(colData(spe)$SignalDensity),
    target_sum = as.numeric(colData(spe)$target_sum),
    control_sum = as.numeric(colData(spe)$control_sum),
    total = as.numeric(colData(spe)$total),
    formula_variables = metadata(spe)$formula_variables,
    n_cells = ncol(spe),
    n_genes = nrow(spe)
)

# Write JSON
jsonlite::write_json(result, output_path, auto_unbox=TRUE, digits=NA, pretty=TRUE)

cat("R reference completed. Output:", output_path, "\n")
cat("QC_score range:", range(result$QC_score), "\n")
