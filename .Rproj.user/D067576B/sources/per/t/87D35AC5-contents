data(fung)

# prep data
fung_detect <- fung %>%
  dplyr::select(1:4)

fung_sitecov_tbl <- fung %>%
  dplyr::select(-sample, -pcr1, -pcr2) %>%
  unique(.)

fung_sampcov_tbl <- fung %>%
  dplyr::select(-pcr1, -pcr2)

# fit intercept model at all three levels use beta-binomial sampler
fung_mod1 <- MSOcc_mod(wide_data = fung_detect, progress = T,
                       site = list(model = ~ 1, cov_tbl = fung_sitecov_tbl),
                       sample = list(model = ~ 1, cov_tbl = fung_sampcov_tbl),
                       num.mcmc = 5000, beta_bin = T)

theta_mcmc(fung_mod1)
