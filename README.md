# Monte Carlo Fair Value Analysis for Loan Pools — Companion Code

Companion code and sample data for the book *Monte Carlo Fair Value Analysis
for Loan Pools: Dynamic Nelson-Siegel Yield Curves, Stochastic Simulation,
Prepayment Modeling, and OAS Analytics*, by Ariel Rodríguez Herrera
(CEO, Ambient Point Corp.), first edition, 2026.

ISBN 979-8-234-25206-7 (print) · ISBN 979-8-234-25207-4 (ebook)

## About the book

The book develops the technical methodology behind Monte Carlo fair value
analysis for pools of amortizing loans and their associated interest-only (IO)
and servicing strips. The treatment is self-contained: it begins with yield
curve fundamentals, builds up the Dynamic Nelson-Siegel (DNS) term structure
model, explains how to simulate interest rate scenarios via stochastic
differential equations, develops at length the theory of prepayment and credit
behavior that drives the economics of these assets, and then applies the
simulated paths to project loan cashflows and compute option-adjusted fair
value.

The emphasis is on concepts and mathematics first, with working code
illustrations second. The methodology is general and applies to any Monte Carlo
valuation engine for amortizing assets; no particular software product is
required to follow the exposition. The reader is assumed to have working
familiarity with undergraduate probability and calculus and with basic
fixed-income concepts such as yield, discount factor, and duration.

**Part I — The Term Structure and the Nelson-Siegel Model** (Chapters 1–6)
lays the mathematical and financial foundations: yield curve theory, the
Nelson-Siegel decomposition, its dynamic extension, stochastic differential
equations, numerical simulation schemes, and the YUIMA framework.

**Part II — Prepayment and Credit Modeling** (Chapters 7–13) develops the
theory and practice of prepayment and credit modeling — the analytical core of
the book — through survival analysis, econometric models, machine-learning
approaches, a loan pool credit model, and validation.

**Part III — Cashflow Projection** (Chapters 14–18) covers amortization
mechanics, the cashflow components, settlement and accrued interest, pool
aggregation, and diagnostics.

**Part IV — Monte Carlo Valuation and Risk Analytics** (Chapters 19–22)
presents path-level NPV and convergence diagnostics, IO and servicing strip
valuation, OAS analysis, sensitivity and key rate duration, and hedging with
Treasury futures.

The **appendices** give an overview of the Q-ABS Fair Value Engine — the
production implementation of this methodology — together with a function
reference and a glossary.

## What is in this repository

```
mc-fair-value-book/
├── 01-term-structure.R  …  22-hedging.R   One file per chapter
├── data/
│   ├── README.md                          Provenance and column reference
│   ├── make_sample_data.R                 Regenerates the two panels
│   ├── sample_loan_tape.csv               135 loans (pool extract)
│   ├── sample_loan_panel.csv              1,500 loans (survival format)
│   └── sample_monthly_panel.csv           53,550 loan-months
└── LICENSE
```

Each `.R` file holds every R code block from the corresponding chapter, in
order, with a banner naming the section and the line of the manuscript source
it came from:

```r
## ---- 01: CPR ⇄ SMM Conversions  (qmd line 238) ----
```

These files are generated from the book manuscript. Corrections belong in the
book source rather than here; edits made directly to these files are lost on
the next regeneration.

## Chapter map

### Part I — The Term Structure and the Nelson-Siegel Model

| File | Chapter |
| :--- | :--- |
| `01-term-structure.R` | The Term Structure Problem |
| `02-nelson-siegel.R` | The Nelson-Siegel Model |
| `03-dynamic-ns.R` | The Dynamic Nelson-Siegel Model |
| `04-sde-simulation.R` | Stochastic Differential Equations and Numerical Simulation |
| `05-yuima.R` | The YUIMA Framework and DNS Simulation |
| `06-yield-evolution.R` | Yield Curve Evolution: From Factors to Scenarios |

### Part II — Prepayment and Credit Modeling

| File | Chapter | Sample data |
| :--- | :--- | :--- |
| `07-prepayment-foundations.R` | Prepayment and Default Modeling: Foundations | loan panel |
| `08-prepayment-theory.R` | The Theory of Prepayment Behavior | — |
| `09-survival-analysis.R` | Survival Analysis and Hazard Models | loan panel, monthly panel |
| `10-econometric-models.R` | Econometric Prepayment Models | monthly panel |
| `11-ml-prepayment.R` | Machine Learning Approaches to Prepayment | monthly panel |
| `12-credit-model.R` | The Loan Pool Credit Model | loan panel |
| `13-model-validation.R` | Model Validation and Benchmarking | loan panel, monthly panel |

### Part III — Cashflow Projection

| File | Chapter | Sample data |
| :--- | :--- | :--- |
| `14-cashflow-projection.R` | Cashflow Projection | — |
| `15-amortization-mechanics.R` | Amortization Mechanics | — |
| `16-settlement-conventions.R` | Settlement Conventions | — |
| `17-pool-aggregation.R` | Pool Aggregation | loan tape |
| `18-cashflow-diagnostics.R` | Cashflow Diagnostics | — |

### Part IV — Monte Carlo Valuation and Risk Analytics

| File | Chapter |
| :--- | :--- |
| `19-fair-value-mc.R` | Monte Carlo Fair Value |
| `20-oas-analysis.R` | Option-Adjusted Spread Analysis |
| `21-sensitivity.R` | Sensitivity Analytics |
| `22-hedging.R` | Hedging |

## Running the code

R 4.0 or later. Clone the repository and run from its root, so that the
relative data paths in the chapter code (`data/sample_loan_panel.csv`, and so
on) resolve.

Within a chapter, blocks are written to run in order: later blocks may use
objects defined by earlier ones. The Part IV chapters (19 through 22) each
open with a shared setup block that builds the pool, rate paths, bond-level
cashflows and solved OAS from the sample pool bundled with `oasr`, so each
chapter runs on its own. Chapter 20's reporting block calls
`create_oas_report()`, which renders through pandoc: run it from RStudio, or
point `RSTUDIO_PANDOC` at a pandoc installation first.

### Packages

The examples run on CRAN packages together with two MIT-licensed packages,
`dnsr` (DNS calibration and rate-path simulation) and `oasr` (cashflow
projection and OAS analytics):

```r
# install.packages("remotes")
remotes::install_github("ambientpointcorp/dnsr")
remotes::install_github("ambientpointcorp/oasr")
```

From CRAN, by part:

```r
# Part I — term structure and simulation
install.packages(c("dplyr", "tidyr", "ggplot2", "KFAS", "vars", "yuima"))

# Part II — prepayment and credit modeling
install.packages(c("survival", "survminer", "cmprsk", "flexsurv", "survAUC",
                   "actuar", "MASS", "broom", "glmnet", "gratia", "mgcv",
                   "mfx", "nnet", "plm", "sandwich", "lmtest", "pROC"))

# Chapter 11 — machine learning (heavier; skip if not needed)
install.packages(c("tidymodels", "xgboost", "lightgbm", "bonsai", "ranger",
                   "shapviz"))

# Parts III and IV — cashflows, valuation, risk
install.packages(c("bizdays", "future.apply", "randtoolbox", "profvis",
                   "PerformanceAnalytics", "digest", "purrr", "readr",
                   "tibble", "zoo"))
```

### FRED API key

Part I pulls a Treasury yield panel from FRED through `dnsr`. A free API key is
available at <https://fred.stlouisfed.org/docs/api/api_key.html>; put it in the
`fred_api_key` environment variable — for example in `~/.Renviron`:

```
fred_api_key=your_key_here
```

### Blocks that need your own inputs

Parts I and II run end to end against the bundled sample data. In Parts III and
IV, the blocks that call the valuation pipeline operate on objects built
upstream in the book's workflow — an aggregated cashflow list from
`aggregate_cashflows()`, a rate-paths object from `simulate_rates()` — rather
than on bundled data. Each such block names its inputs in a comment at the top:

```r
# pool_cf: list returned by aggregate_cashflows()
# rate_paths: result of simulate_rates(); $yield_paths holds
#   [n_paths x n_periods] matrices keyed by maturity
```

Substitute the equivalent objects for your own pool. The chapters show how each
is constructed.

## Sample data

Three files support the empirical examples; `data/README.md` documents them in
full, with a column reference for each.

`sample_loan_tape.csv` is an anonymized extract from a securitized pool of
Prime-indexed, government-guaranteed small-business loans. Row-level
characteristics — terms, spreads, the fee stack, balances, sectors, states —
are genuine; identifying attributes were removed, dates and balances jittered,
NAICS codes coarsened to two-digit sectors, and row order shuffled.

`sample_loan_panel.csv` and `sample_monthly_panel.csv` are simulated by
`data/make_sample_data.R`: loan characteristics are bootstrapped from the tape,
and monthly prepayment and default histories are drawn from seasoning,
refinance-incentive, sector, and vintage hazards calibrated to the aggregate
behavior of the public SBA 7(a) FOIA dataset. The histories are synthetic; no
borrower-level performance data is included. To regenerate them:

```bash
cd data && Rscript make_sample_data.R
```

## License

The code and sample data in this repository are released under the MIT License
(see `LICENSE`). You may use these examples in your own programs and
documentation, including commercially, without seeking permission. Permission
*is* required to reproduce or redistribute a significant portion of them as
such — for instance publishing a collection of this code, or incorporating
substantial portions of it into another book, course, or commercial training
product.

Attribution is appreciated but not required: *Monte Carlo Fair Value Analysis
for Loan Pools*, © 2026 Ariel Rodríguez Herrera, ISBN 979-8-234-25207-4.

The book's text, figures, and tables are © 2026 Ariel Rodríguez Herrera, all
rights reserved. The `dnsr` and `oasr` packages are open source under the MIT
License. Ambient Point Corp.'s commercial software — including the Q-ABS Fair
Value Engine and the Q-SBA suite of packages and applications — is proprietary
and licensed separately.

The code, models, and procedures here are provided "as is," without warranty of
any kind, for instruction rather than as investment, accounting, or legal
advice. Nothing here is a recommendation to buy or sell any security.

## Links

- [ambientpointcorp.com](https://ambientpointcorp.com)
- [quantome.ai](https://quantome.ai)
