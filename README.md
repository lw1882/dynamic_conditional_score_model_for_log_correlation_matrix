# dynamic conditional score model for log correlation matrix

Code to estimate the dynamic conditional score (DCS) model for the **log correlation matrix** proposed in:

- Hafner, C. M., and Wang, L. (2023). *A dynamic conditional score model for the log correlation matrix*.

**Input:** Sample daily return data stored in `data/asset_returns.csv`.

**Model estimation**
- Step 1: Run `main_01_volatility.R` to estimate volatilities using the Beta-t-EGARCH model.
- Step 2: Run `main_02_correlation.R` to estimate correlations using the log-correlation DCS model.
