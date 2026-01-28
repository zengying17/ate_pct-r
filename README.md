# ate_pct (R)

R implementation of **ATE in percentage points (ATE-pct)** under subgroup heterogeneity, as developed in:

> Ying Zeng, *“Estimation and Inference on Average Treatment Effect in Percentage Points”*.

This repository provides:
- `ate_pct()` — post-estimation function for `lm()`-like objects (requires `coef()` and `vcov()` methods)
- an example dataset (`ate_pct_example.Rda`)

## Contents

- `ate_pct.R` — main function and methods
- `ate_pct_example.Rda` — example dataset

## Installation

### Option A: Source directly from GitHub (download + source)
1. Download or clone this repository.
2. In R:
   ```r
   source("ate_pct.R")
   ```

### Option B: Use `source()` from a local path
```r
source("/path/to/ate_pct-r/ate_pct.R")
```

## Quick start

```r
library("sandwich")

# Load functions
source("ate_pct.R")

# Load example data
load("ate_pct_example.Rda")

# Example regression (adjust variable names to your dataset)
reg_res <- lm(lny~x+gr1+gr2+gr3,data=df)
VChet <- vcovHC(reg_res, type = "HC1")
# Post-estimation
out <- ate_pct(reg_res, c("gr1","gr2","gr3"), VC_model=VChet)
print(out)
summary(out)
```
See `ate_pct.R` for additional worked examples and usage notes (included in the script comments).

## Notes on sampling and variance output

- Subsetting should be handled when fitting the model on the desired sample (e.g., via `subset=` in `lm()`), because `ate_pct()` uses the fitted model’s coefficients and `vcov()`.
- `vcov(out)` returns a 3×3 **diagonal** matrix containing delta-method variances for `(taubar, rho_a, rho_b)`; off-diagonal covariances are ignored and set to 0 by design.

## Citation

If you use this code in academic work, please cite:

- Zeng, Ying. *Estimation and Inference on Average Treatment Effect in Percentage Points.*

## Contact

Ying Zeng  

School of Economics and Wang Yanan Institute for Studies in Economics (WISE)

Xiamen University

Email: <zengying17@gmail.com>

## Related repository

Stata implementation: [ate_pct-stata](https://github.com/zengying17/ate_pct-stata)


