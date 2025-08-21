
test_that("fair_surv: copula sensitivity mode yields tau-stratified SURV measures", {
  dt <- mk_data(100L)  # keep competing-risk status
  X <- "majority"; Z <- "age"; W <- c("clinstg","hgb","ch")
  
  taus <- c(0.3, 0.8)
  fs <- fair_surv(
    data = dt, X = X, Z = Z, W = W,
    time_var = "time", event_var = "status",
    nboot = 3,
    copula = "clayton",
    tau_grid = taus
  )
  
  expect_false(fs$is_cr)   # in sensitivity, computation path forces nslices = 1
  expect_true(fs$is_sens)
  
  # measures restricted to 'surv' and include tau/copula columns
  expect_setequal(unique(fs$measures$scale), "surv")
  expect_true(all(c("tau","copula") %in% names(fs$measures)))
  expect_setequal(unique(fs$measures$tau), taus)
  expect_setequal(unique(fs$measures$copula), "clayton")
  
  # event column should be 1 (no slicing by cause in sensitivity loop)
  expect_setequal(unique(fs$measures$event), 1L)
  
  # values finite and sd present from bootstrap aggregation
  expect_true(all(is.finite(fs$measures$value)))
  expect_true(all(is.finite(fs$measures$sd)))

  for (sc in c("surv")) {
    p <- autoplot(fs, scale=sc)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot_build(p))
  }
})
