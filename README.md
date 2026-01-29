# ate_pct (R)

R implementation of **ATE in percentage points (ATE-pct)** under subgroup heterogeneity, as developed in:

> Ying Zeng, *“Estimation and Inference on Average Treatment Effect in Percentage Points under Heterogeneity”*.

This repository provides:
- `ate_pct()` — post-estimation function for `lm()`-like objects (requires `coef()` and `vcov()` methods)
- an example dataset (`ate_pct_example.Rda`)

## Contents

- `ate_pct.R` — main function and methods
- `ate_pct_example.Rda` — example dataset

## Installation
```
install.packages("remotes")
remotes::install_github("zengying17/ate_pct-r")
```


## Quick start

```r
library("sandwich")
library("atepct")
# Load example data

data(ate_pct_example, package = "atepct")

# Example regression (adjust variable names to your dataset)
reg_res <- lm(lny~x+gr1+gr2+gr3,data=ate_pct_example)
# Default: uses vcov(reg_res)
out1 <- ate_pct(reg_res, c("gr1","gr2","gr3"))
out1
summary(out1)

# Optional: heteroskedasticity-robust vcov via sandwich
if (requireNamespace("sandwich", quietly = TRUE)) {
  VChet <- sandwich::vcovHC(reg_res, type = "HC1")
  out1_robust <- ate_pct(reg_res, c("gr1","gr2","gr3"), VC_model = VChet)
  summary(out1_robust)
}
```
See the help document for additional worked examples and usage notes.

## Notes on sampling and variance output

- Subsetting should be handled when fitting the model on the desired sample (e.g., via `subset=` in `lm()`), because `ate_pct()` uses the fitted model’s coefficients and `vcov()`.
- `vcov(out)` returns a 3×3 **diagonal** matrix containing delta-method variances for `(taubar, rho_a, rho_b)`; off-diagonal covariances are ignored and set to 0 by design.

## Included example data

The package includes ate_pct_example for replicable examples:

```r
data(ate_pct_example, package = "atepct")
```

## Citation

If you use this code in academic work, please cite:

- Zeng, Ying. *Estimation and Inference on Average Treatment Effect in Percentage Points under Heterogeneity.*

## Contact

Ying Zeng  

School of Economics and Wang Yanan Institute for Studies in Economics (WISE)

Xiamen University

Email: <zengying17@gmail.com>

## Related repository

Stata implementation: [ate_pct-stata](https://github.com/zengying17/ate_pct-stata)


