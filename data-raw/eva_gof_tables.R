# One-off script to pull internal GOF quantile tables from eva package
# and store them as internal data in permApprox.

library(eva)

# Grab unexported objects from eva's namespace
ADQuantiles_eva  <- get("ADQuantiles",  envir = asNamespace("eva"))
CVMQuantiles_eva <- get("CVMQuantiles", envir = asNamespace("eva"))

# Give them internal names for permApprox
.ADQ_gpd  <- ADQuantiles_eva
.CVMQ_gpd <- CVMQuantiles_eva

# Save as internal data (R/sysdata.rda)
usethis::use_data(
  .ADQ_gpd,
  .CVMQ_gpd,
  internal  = TRUE,
  overwrite = TRUE
)
