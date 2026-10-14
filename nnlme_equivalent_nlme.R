# nnlme_equivalent_nlme.R
# R translation of the Stata command nnlme (Sean Clouston; SSC S459824, revised 10 Oct 2026; in Stata: ssc install nnlme).
#
# nnlme fits, with Stata's menl:
#   y_ij  = {ri1:} + {shp:} + (-exp(a00)) * ln(1 + exp(min(700, time_ij - exp(-eta_i))))
#   ri1   = gamma_0 + U0_i                (intercept plus random intercept: define(ri1: U0[id], xb))
#   shp   = gamma' w_ij                   (shape() variables, no constant; without shape() this term is absent)
#   eta_i = beta_0 + beta' x_i            (surv() variables, with a constant); onset time tau_i = exp(-eta_i)
#   U0_i ~ N(0, sd_u^2), e_ij ~ N(0, sigma^2)
# Here the same model is fitted with nlme::nlme(), which, like menl, uses the Lindstrom-Bates linearization and maximum likelihood.
# nlme attaches a random effect to a fixed parameter, so gamma_0 and U0_i sit on the constant of shp; this is the same model as {ri1:} + {shp:} in the Stata command.
# nnlme shows exp(beta) for the surv() variables as NIR: one unit higher x multiplies the onset time by 1/NIR.
#
# Usage: edit the settings in section 0, then run with the working directory set to this folder: source("nnlme_equivalent_nlme.R")
# Requires tidyverse and nlme. On the simulated data it runs in well under a minute.
# Stata original: ssc install nnlme (https://ideas.repec.org/c/boc/bocode/s459824.html). On the same data the two should agree closely; small differences can come from the optimizers.

suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
})

# 0. Settings: same names as the nnlme options ------------------
# The values below reproduce  nnlme y, surv(groupb groupc z) id(id) time(time) shape(time x z groupb groupc)
data_file <- "simulated_onset_data.csv"               # CSV in long format, one row per person per visit
depvar    <- "y"                                       # outcome (nnlme depvar)
id        <- "id"                                      # person identifier, id()
time      <- "time"                                    # time variable, time(); its median must be > 0 (used for the starting value)
surv      <- c("groupb", "groupc", "z")                # surv(): variables that shift the onset time, shown as NIR
shape     <- c("time", "x", "z", "groupb", "groupc")   # shape(); character(0) for nnlme without shape()
iter      <- 25                                        # iter(); nnlme's default is 25

# 1. Data ------------------
# Rows with a missing value in any variable of the model are dropped, as Stata does
d_raw <- read_csv(data_file, show_col_types = FALSE)
d <- d_raw %>%
  select(all_of(unique(c(id, depvar, time, surv, shape)))) %>%
  drop_na()
cat(sprintf("rows %d (%d dropped for missing values), people %d; %s %.1f to %.1f\n",
            nrow(d), nrow(d_raw) - nrow(d), n_distinct(d[[id]]), time, min(d[[time]]), max(d[[time]])))

# 2. Model ------------------
# Formulas are built from the settings; factor or character columns get dummies with the first level as the base, like Stata's i. with the lowest value
# Starting values come from the data only: intercept = mean outcome, onset time = median time, decline 0.5 per year, all other coefficients 0
f_model <- as.formula(sprintf("%s ~ shp - exp(a00) * log(1 + exp(pmin(700, %s - exp(-eta))))", depvar, time))
f_shp <- reformulate(if (length(shape) > 0) shape else "1", response = "shp")
f_eta <- reformulate(surv, response = "eta")
f_random <- as.formula(sprintf("shp ~ 1 | %s", id))
n_shp <- ncol(model.matrix(delete.response(terms(f_shp)), d))
n_eta <- ncol(model.matrix(delete.response(terms(f_eta)), d))

fit <- nlme(f_model,
            data = as.data.frame(d),
            fixed = list(f_shp, a00 ~ 1, f_eta),
            random = f_random,
            start = c(shp = c(mean(d[[depvar]]), rep(0, n_shp - 1)), a00 = log(0.5), eta = c(-log(median(d[[time]])), rep(0, n_eta - 1))),
            method = "ML",
            control = nlmeControl(maxIter = iter))

# 3. Results in nnlme's layout ------------------
# Blocks, row names and order as in the SSC nnlme output. For the /survival rows nnlme shows NIR = exp(coef), Std. Err. = NIR * SE (delta method)
# and the CI as exp() of the log-scale CI; z and p use the normal distribution on the log scale
tt <- summary(fit)$tTable
res <- tibble(nlme_term = rownames(tt), coef = tt[, "Value"], se = tt[, "Std.Error"]) %>%
  mutate(block = case_when(str_starts(nlme_term, "eta\\.") ~ "/survival",
                           nlme_term == "shp.(Intercept)" ~ "Intercept",
                           str_starts(nlme_term, "shp\\.") ~ "/shape",
                           nlme_term == "a00" ~ "Acc. param."),
         block = factor(block, levels = c("/survival", "Intercept", "/shape", "Acc. param.")),
         term = nlme_term %>% str_remove("^(eta|shp)\\.") %>% str_replace("^\\(Intercept\\)$", "_cons"),
         is_surv = block == "/survival",
         z_stat = coef / se,
         p = 2 * pnorm(-abs(z_stat)),
         shown = if_else(is_surv, exp(coef), coef),
         shown_se = if_else(is_surv, exp(coef) * se, se),
         ci_low = if_else(is_surv, exp(coef - 1.96 * se), coef - 1.96 * se),
         ci_high = if_else(is_surv, exp(coef + 1.96 * se), coef + 1.96 * se)) %>%
  arrange(block, term == "_cons")
sd_u <- as.numeric(VarCorr(fit)[1, "StdDev"])

cat(sprintf("\nNested Non-Linear Longitudinal Mixed-Effects Model (R, nlme)\nNumber of obs %d; groups (%s) %d; log likelihood %.2f\n\n",
            nobs(fit), id, n_distinct(d[[id]]), as.numeric(logLik(fit))))
res %>%
  transmute(block = as.character(block),
            term,
            `NIR / Coef.` = sprintf("%.7g", shown),
            `Std. Err.` = sprintf("%.7g", shown_se),
            z = sprintf("%.2f", z_stat),
            `P>|z|` = sprintf("%.3f", p),
            `95% CI low` = sprintf("%.7g", ci_low),
            `95% CI high` = sprintf("%.7g", ci_high)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
cat(sprintf("\nRandom Eff.\n/%s lnsd(U0) %.7g (SD %.7g)\n/residual lnsigma %.7g (SD %.7g)\n",
            id, log(sd_u), sd_u, log(fit$sigma), fit$sigma))

# 4. Simulated data only: compare with the true values used in simulate_onset_data.R ------------------
if (basename(data_file) == "simulated_onset_data.csv") {
  truth <- c(`shp.(Intercept)` = 25, shp.time = -0.05, shp.x = 1, shp.z = 0.3, shp.groupb = 0, shp.groupc = 0,
             a00 = log(0.6),
             `eta.(Intercept)` = -log(15), eta.groupb = -log(22 / 15), eta.groupc = 0, eta.z = -log(16 / 15))
  eta_hat <- unname(fixef(fit)[c("eta.(Intercept)", "eta.groupb", "eta.groupc")])
  onset_age <- 50 + exp(-(eta_hat[1] + c(A = 0, B = eta_hat[2], C = eta_hat[3])))
  cat("\nSimulated truth (truth_shown on nnlme's display scale; dist_in_se = (estimate - truth) / SE on the estimation scale):\n")
  res %>%
    transmute(block = as.character(block),
              term,
              shown = sprintf("%.4f", shown),
              truth_shown = sprintf("%.4f", if_else(is_surv, exp(truth[nlme_term]), truth[nlme_term])),
              dist_in_se = sprintf("%.2f", (coef - truth[nlme_term]) / se)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
  cat(sprintf("Random intercept SD %.4f (truth 2); residual SD %.4f (truth 1.5)\n", sd_u, fit$sigma))
  cat(sprintf("Onset age at z = 0 (50 + exp(-eta)): A %.1f, B %.1f, C %.1f (truth 65, 72, 65)\n", onset_age["A"], onset_age["B"], onset_age["C"]))
  cat(sprintf("Decline after the onset: %.3f per year (truth -0.600)\n", -exp(fixef(fit)["a00"])))
}
cat(sprintf("\n%s; nlme %s\n", R.version.string, as.character(packageVersion("nlme"))))
