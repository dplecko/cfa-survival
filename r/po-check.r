# 
# poutx1 <- fsurv$srv$srvx1[, 100]
# poutx0 <- fsurv$srv$srvx0[, 100]
# pout <- fsurv$srv$srvx[, 100]
# 
# px_zw <- fsurv$pw$px_zw
# px_z <- fsurv$pw$px_z
# px <- fsurv$pw$px
# 
# X <- "majority"
# 
# compute_po <- function(fx = 0, wx = 0, zx = 0, thresh = FALSE, 
#                        norm_method = "adapt") {
#   
#   if (is.na(zx)) zx <- -1
#   
#   if (wx == 0) wgh <- (1 - px_zw) / (1 - px_z) else wgh <- px_zw / px_z
#   if (zx == 0) {
#     wgh <- wgh * (1 - px_zw) / (1 - px)
#   } else if (zx == 1) wgh <- wgh * px_z / px
#   
#   if (fx == 0) po_samp <- poutx0 else po_samp <- poutx1
#   if (norm_method == "adapt") norm_const <- sum(wgh) else 
#     norm_const <- length(po_samp)
#   if (thresh) po_samp <- as.integer(po_samp > thr)
#   list(po_samp = po_samp, wgh = wgh, norm_const = norm_const)
# }
# 
# eval_po <- function(po) sum(po$po_samp * po$wgh) / po$norm_const
# diff_po <- function(po1, po2) eval_po(po1) - eval_po(po2)
# 
# norm_method = "adapt"
# 
# f_x1 <- mean(pout[dat_run[[X]] == 1])
# 
# fx1wx1_x1 <- compute_po(fx = 1, wx = 1, zx = 1, FALSE, norm_method)
# eval_po(fx1wx1_x1)
# 
# f_x0 <- mean(pout[dat_run[[X]] == 0])
# fx0wx0_x0 <- compute_po(fx = 0, wx = 0, zx = 0, FALSE, norm_method)
# 
# eval_po(fx0wx0_x0)
