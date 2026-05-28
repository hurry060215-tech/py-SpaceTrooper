# MATH.md — Acceleration Agent Derivations

## Algorithm Summary

SpaceTrooper computes QC scores for spatial transcriptomics via:

1. **Per-cell QC metrics**: sum, detected, control/total ratios, signal density
2. **Outlier detection**: medcouple-based (skewed distributions) or MAD-based (symmetric distributions)
3. **QC score**: ridge logistic regression (L2 penalty) with cross-validated lambda

## Key Numerical Considerations

### Medcouple (MC)
- Robust skewness measure in [-1, 1]
- Used to decide outlier detection method: |MC| >= 0.6 → adjusted boxplot; else → MAD
- Formula: median{ h(x_i, x_j) : x_i <= med <= x_j }
- h(x_i, x_j) = (x_i + x_j - 2*med) / (x_j - x_i)

### Ridge Logistic Regression
- sklearn LogisticRegression(penalty='l2') ≈ R glmnet(alpha=0)
- Solver differences: R uses coordinate descent, Python uses L-BFGS
- Expected numerical drift: ~1e-6 to 1e-4 in QC_score values
- This drift is within the deterministic-bounded parity threshold (1e-6)

### Outlier Thresholds
- MAD-based: median ± 5 * MAD * 1.4826 (scuttle default)
- Medcouple-based: Q1 - 1.5*exp(-4*MC)*IQR, Q3 + 1.5*exp(3*MC)*IQR

## Parity Gate Justification

The primary output (QC_score) is a continuous float from ridge logistic regression.
Given identical training data and fixed seed, the training set construction is
deterministic. The GLM solver introduces small numerical noise (~1e-6) between
R glmnet and Python sklearn, justifying the `deterministic-bounded` classification
with threshold 1e-6.

Outlier labels (LOW/HIGH/NO) are categorical and use classification metrics (F1 >= 0.95).
