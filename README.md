# TGP_for_g-2

_Fit the R-ratio using the Treed Guassian process（TGP） model and perform a series of calculations._

## Fitting:

_The TGP model was applied to the collated dataset, and the GP model was fitted to enable comparisons._

```bash
TGP/TGP_for_R_ratio.py
TGP/GP_for_R_ratio.py
```

## Calculation with TGP:

_Calculate the values of amuon and alpha, along with their uncertainties, and explore their relationship using TGP. Additionally, search for the location of maximum uncertainty._

```bash
TGP/Compute_with_TGP.py
TGP/Cov_with_TGP.py
TGP/ALM_for_amuon.py
TGP/ALM_for_alpha.py
```

## Calculation with naive model:

_Calculate the values and uncertainties of amuon and alpha solely based on data._

```bash
TGP/Naive_calculation_for_amuon.py
TGP/Naive_calculation_for_alpha.py
TGP/Naive_calculation_for_amuon_uncertainty.py
TGP/Naive_calculation_for_alpha_uncertainty.py
```

&starf;Feel free to explore the respective files for detailed implementation and usage instructions.
