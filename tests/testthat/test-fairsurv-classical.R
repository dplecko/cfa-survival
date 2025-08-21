
test_that("fair_surv: classical survival mode works and structure is valid", {
  dt <- mk_data(100L)
  # collapse to 1 event
  dt[, status := as.integer(status > 0)]
  
  X <- "majority"; Z <- "age"; W <- c("clinstg","hgb","ch")
  fs <- fair_surv(
    data = dt, X = X, Z = Z, W = W,
    time_var = "time", event_var = "status",
    nboot = 3
  )
  
  expect_s3_class(fs, "fairsurv")
  expect_false(fs$is_cr)
  expect_false(fs$is_sens)
  
  # measures table basics
  expect_true(is.data.table(fs$measures))
  expect_setequal(
    names(fs$measures),
    c("time_interest","effect","scale","event","value","sd")
  )
  expect_setequal(unique(fs$measures$scale),
                  c("surv","mst","chf","chf-ratio"))
  expect_equal(length(unique(fs$measures$event)), 1L)
  expect_true(all(is.finite(fs$measures$value)))
  expect_true(all(is.finite(fs$measures$sd)))
  
  # time grid wired through
  expect_true(all(fs$time_interest %in% fs$measures$time_interest))
  
  # survival trajectories should be non-increasing over time (few rows check)
  s <- fs$srv$srvx
  expect_true(ncol(s) == length(fs$time_interest))
  for (i in 1:5) {
    expect_true(all(diff(s[i, ]) <= 1e-8))
    expect_true(all(s[i, ] >= 0 & s[i, ] <= 1))
  }
  
  # CIFs absent in classical mode
  expect_null(fs$cif$cifx)
  expect_null(fs$cif$cif_x0)
  expect_null(fs$cif$cif_x1)
  
  # CHF dims coherent
  chf <- fs$chf$chfx
  expect_equal(dim(chf), c(nrow(dt), length(fs$time_interest)))
  
  for (sc in c("surv", "mst", "chf", "chf-ratio")) {
    p <- autoplot(fs, scale=sc)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot_build(p))
  }
})
