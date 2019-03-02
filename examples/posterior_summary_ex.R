data(fung)

# prep data
fung_detect <- fung %>%
  dplyr::select(1:4)

fung_sitecov_tbl <- fung %>%
  dplyr::select(-sample, -pcr1, -pcr2) %>%
  unique(.)

fung_sampcov_tbl <- fung %>%
  dplyr::select(-pcr1, -pcr2)

# model sample level occurence by frog density
fung_mod2 <- MSOcc_mod(wide_data = fung_detect, progress = T,
                       site = list(model = ~ 1, cov_tbl = fung_sitecov_tbl),
                       sample = list(model = ~ frogs, cov_tbl = fung_sampcov_tbl),
                       num.mcmc = 5000, beta_bin = T)

# overall summary
posterior_summary(fung_mod2)

# print unique rows
posterior_summary(fung_mod2) %>%
  dplyr::mutate(site = unlist(lapply(strsplit(site_sample, '_'), function(x) x[1]))) %>%
  dplyr::select(-site_sample) %>%
  dplyr::select(site, dplyr::everything()) %>%
  dplyr::distinct()

# site level
posterior_summary(fung_mod2, level = 'site')

# sample level
posterior_summary(fung_mod2, level = 'sample')

# rep level
posterior_summary(fung_mod2, level = 'rep')

