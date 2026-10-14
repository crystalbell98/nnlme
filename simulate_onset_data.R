# simulate_onset_data.R
# Builds simulated_onset_data.csv, the demo data used by nnlme_equivalent_nlme.R.
# The CSV is shipped with the scripts, so re-running this file is optional; it documents how the data were made and the true values.
#
# Data-generating model = the nnlme model with shape(), with time = age - 50:
#   y_ij  = U0_i + shp_ij + (-exp(a00)) * ln(1 + exp(time_ij - tau_i)) + e_ij
#   shp   = 25 - 0.05 * time + 1 * x + 0.3 * z + 0 * groupb + 0 * groupc
#   tau_i = exp(-eta_i), eta_i = -ln(15) - ln(22 / 15) * groupb + 0 * groupc - ln(16 / 15) * z
#   a00 = ln(0.6), i.e. 0.6 points per year faster decline after the onset
#   U0_i ~ N(0, 2^2), e_ij ~ N(0, 1.5^2)
# What these values mean:
#   onset age at z = 0: group A 50 + 15 = 65, group B 50 + 22 = 72, group C 65 (same as A on purpose)
#   NIR (exp of the eta coefficients): groupb 15/22 = 0.682, groupc 1, z 15/16 = 0.9375 (each unit of z multiplies the onset time by 16/15)
# Design: 400 people, 5 tests each, 2 years apart, age at entry uniform 50-82; x = 0/1 (e.g., sex), z ~ N(0, 3^2) (e.g., centred years of education).

library(tidyverse)

set.seed(20261013)
N_PEOPLE <- 400
people <- tibble(id = seq_len(N_PEOPLE),
                 group = sample(c("A", "B", "C"), N_PEOPLE, replace = TRUE),
                 x = rbinom(N_PEOPLE, 1, 0.5),
                 z = rnorm(N_PEOPLE, 0, 3),
                 age_entry = runif(N_PEOPLE, 50, 82),
                 u0 = rnorm(N_PEOPLE, 0, 2))
d <- people %>%
  crossing(wave = 0:4) %>%
  mutate(groupb = as.integer(group == "B"),
         groupc = as.integer(group == "C"),
         age = age_entry + 2 * wave,
         time = age - 50,
         eta = -log(15) - log(22 / 15) * groupb + 0 * groupc - log(16 / 15) * z,
         y = 25 - 0.05 * time + 1 * x + 0.3 * z + u0 - exp(log(0.6)) * log(1 + exp(time - exp(-eta))) + rnorm(n(), 0, 1.5)) %>%
  select(id, group, groupb, groupc, x, z, age, time, y)

write_csv(d, "simulated_onset_data.csv")
cat(sprintf("wrote simulated_onset_data.csv: %d rows, %d people\n", nrow(d), n_distinct(d$id)))
