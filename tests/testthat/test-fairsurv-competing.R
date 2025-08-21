
test_that("fair_surv: competing risks mode exposes CIF and multi-event measures", {
  dt <- mk_data(100L) # original status has >1 event
  X <- "majority"; Z <- "age"; W <- c("clinstg","hgb","ch")
  
  fs <- fair_surv(
    data = dt, X = X, Z = Z, W = W,
    time_var = "time", event_var = "status",
    nboot = 3
  )
  
  expect_true(fs$is_cr)
  expect_false(fs$is_sens)
  
  # scales include CIF
  expect_true("cif" %in% unique(fs$measures$scale))
  
  # event slices cover all causes
  ne <- max(dt$status)
  expect_equal(sort(unique(fs$measures$event)), seq_len(ne))
  
  # CIF arrays present with [n, T, E] shape
  # (cifx is NULL in classical; here should be 3-d arrays)
  # Depending on your implementation, these may be arrays or lists; check robustly.
  # We verify that at least x0/x1 exist and align with time grid.
  expect_true(!is.null(fs$cif$cif_x0))
  expect_true(!is.null(fs$cif$cif_x1))
  
  cx0 <- fs$cif$cif_x0
  cx1 <- fs$cif$cif_x1
  # Accept either 3D array or numeric with dim attr
  dx0 <- dim(cx0); dx1 <- dim(cx1)
  expect_equal(dx0[1], nrow(dt))
  expect_equal(dx0[2], length(fs$time_interest))
  expect_equal(dx0[3], ne)
  expect_equal(dx1, dx0)
  
  # CIFs are in [0,1] and non-decreasing over time for each event
  for (k in seq_len(ne)) {
    v <- cx0[, , k, drop = FALSE]
    # check a few rows
    for (i in 1:5) {
      vi <- as.numeric(v[i, , 1])
      expect_true(all(diff(vi) >= -1e-8))
      expect_true(all(vi >= 0 & vi <= 1))
    }
  }
  
  # expect plotting to work
  for (sc in c("surv", "cif")) {
    p <- autoplot(fs, scale=sc)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot_build(p))
  }
})
