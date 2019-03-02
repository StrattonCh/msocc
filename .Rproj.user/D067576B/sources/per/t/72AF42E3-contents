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

# model sample level occurence by frog density
fung_mod2 <- MSOcc_mod(wide_data = fung_detect, progress = T,
                       site = list(model = ~ 1, cov_tbl = fung_sitecov_tbl),
                       sample = list(model = ~ frogs, cov_tbl = fung_sampcov_tbl),
                       num.mcmc = 5000, beta_bin = T)

# intercept model for site level occurence with random effects for site
fung_mod3 <- MSOcc_mod(wide_data = fung_detect, progress = T,
                       site = list(model = ~ 1 | site, cov_tbl = fung_sitecov_tbl),
                       sample = list(model = ~ 1, cov_tbl = fung_sampcov_tbl),
                       num.mcmc = 5000, beta_bin = T)

# intercept model for site level occurence, model sample level occurence
# by frogs with with random effects for site
fung_mod4 <- MSOcc_mod(wide_data = fung_detect, progress = T,
                       site = list(model = ~ 1 | site, cov_tbl = fung_sitecov_tbl),
                       sample = list(model = ~ frogs | site, cov_tbl = fung_sampcov_tbl),
                       num.mcmc = 5000, beta_bin = T)

