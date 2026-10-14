# nnlme: R module to estimate Nested Nonlinear Longitudinal Mixed-Effects regression

Author: Sean Clouston (Stony Brook University)

An R translation of the Stata command `nnlme`, which estimates the onset of accelerated cognitive decline.
The Stata command is available from SSC ([S459824](https://ideas.repec.org/c/boc/bocode/s459824.html)); install it in Stata with `ssc install nnlme`.
There is no R package; the translation is a single script that uses the CRAN packages nlme and tidyverse.

## Files

| File | What it is |
|---|---|
| `nnlme_equivalent_nlme.R` | Translation of `nnlme` with `nlme::nlme()`; prints results in `nnlme`'s layout, and next to the truth when run on the simulated data |
| `simulated_onset_data.csv` | Demo data: 400 people, 5 tests each (2,000 rows). Columns `id, group, groupb, groupc, x, z, age, time, y`; `time = age - 50` |
| `simulate_onset_data.R` | How the CSV was made, with the true parameter values (re-running it is optional) |

## Usage

Edit the settings in section 0 of the script, then run it with the working directory set to the folder that contains the data:

```r
data_file <- "my_data.csv"
depvar    <- "cogt"
id        <- "idau"
time      <- "cagey"
surv      <- c("raedy", "male")
shape     <- c("stroke", "pstroke", "cagey", "learn")
iter      <- 25

source("nnlme_equivalent_nlme.R")
```

This corresponds to the Stata command:

```stata
nnlme cogt, surv(raedy male) id(idau) time(cagey) shape(stroke pstroke cagey learn)
```

As shipped, the settings run the example on `simulated_onset_data.csv`.

### Arguments

| Setting | nnlme option | Meaning |
|---|---|---|
| `data_file` | (data in memory) | CSV file in long format, one row per person per visit. From Stata: `export delimited using "my_data.csv"` |
| `depvar` | depvar | Outcome variable |
| `id` | `id()` | Person identifier |
| `time` | `time()` | Time variable. Onset time is exp(−η) > 0, so time should be measured from an origin before onset, and its median must be > 0 (it is used for the starting value) |
| `surv` | `surv()` | Variables that shift the onset time. Reported as NIR |
| `shape` | `shape()` | Variables for the mean level and trend. `character(0)` for `nnlme` without `shape()` |
| `iter` | `iter()` | Maximum iterations, passed to `nlmeControl(maxIter =)`. Default 25, as in `nnlme` |

Categorical variables: use 0/1 dummies, or a factor / character column (dummies with the first level as the base, like Stata's `i.` with the lowest value as the base).

### Details

Rows with a missing value in any variable of the model are dropped before fitting, as Stata does.
Starting values come from the data only: intercept = mean outcome, onset time = median time, decline 0.5 per year, all other coefficients 0. Stata uses `menl`'s EM starting values. If the two fits disagree, compare the log likelihoods first: a different value means one of them stopped at a different solution.
`nnlme` passes `pnlsopts(iterate(50) tolerance(1e-8)) emiterate(50)` to `menl`; the R script keeps nlme's own settings for the inner steps and only sets `maxIter`.

### Value (printed output)

| R row | nnlme row |
|---|---|
| `/survival <var>`, `/survival _cons` | `/survival <var>`, `/survival _cons` |
| `Intercept _cons` | `Intercept _cons` |
| `/shape <var>` | `/shape <var>` |
| `Acc. param. a00` | `Acc. param. a00` |
| `/<id> lnsd(U0)` | `/<id> lnsd(U0)` |
| `/residual lnsigma` | `/residual lnsigma` |
| log likelihood | Linearization log likelihood |

As in `nnlme`, the `/survival` rows show NIR = exp(coef), Std. Err. = NIR × SE (delta method) and the CI as exp() of the log-scale CI; z and p use the normal distribution. For the two random-effect rows the R script prints the estimate (and the SD), not its standard error.
When `data_file` is `simulated_onset_data.csv`, the script also prints the true values, onset ages and the decline after the onset.

### Agreement with Stata

Expect coefficients within about 0.1 SE and log likelihoods within about 0.5.
On the simulated data, the SSC `nnlme` (Stata/MP 18.0) and this script gave the same estimates: every coefficient differs by less than 0.01 standard errors, the two SDs by less than 0.001%, and both log likelihoods are -4148.96.
Both recover the simulated truth: onset ages at z = 0 of 65.0, 71.1 and 65.7 years for groups A, B and C (truth 65, 72, 65), NIR 0.707, 0.955 and 0.943 for `groupb`, `groupc` and `z` (truth 0.682, 1 and 0.9375), and a decline of 0.600 per year after the onset (truth 0.600).

### Requirements

R with nlme and tidyverse: `install.packages(c("nlme", "tidyverse"))`. Tested with R 4.5.1 and nlme 3.1.168.

## The model

y_ij = γ0 + U0_i + γ′w_ij − exp(a00) · ln(1 + exp(time_ij − τ_i)) + e_ij, with τ_i = exp(−η_i), η_i = β0 + β′x_i (the `surv()` variables), w_ij the `shape()` variables, U0_i ~ N(0, σ_u²), e_ij ~ N(0, σ²).
NIR = exp(β): one unit higher x multiplies the onset time by 1/NIR. `nlme::nlme()`, like `menl`, fits this by maximum likelihood with the Lindstrom-Bates linearization. In the Stata command the intercept and random intercept are `define(ri1: U0[id], xb)` and the `shape()` terms have no constant; `nlme` attaches the random intercept to the constant of `shp`, which is the same model.

## How to cite

Clouston SAP, Hanes DW, Smith DM. An exploration of sources of bias in population models of cognitive aging and accelerated decline: a prospective cohort study. *J Gerontol B Psychol Sci Soc Sci*. 2026;81(5):gbag040. https://doi.org/10.1093/geronb/gbag040

Sean Clouston, 2026. "NNLME: Stata module to estimate Nested Nonlinear Longitudinal Mixed-Effects regression," Statistical Software Components S459824, Boston College Department of Economics, revised 10 Oct 2026. https://ideas.repec.org/c/boc/bocode/s459824.html

## License

MIT; see [LICENSE](LICENSE).
