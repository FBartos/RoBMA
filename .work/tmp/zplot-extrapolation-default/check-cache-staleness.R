setwd("tests/testthat")
source("common-functions.R")

cat("cache dir:", temp_fits_dir, "\n")
cat("n fit files:", length(list.files(temp_fits_dir, pattern = "[.]RDS$")), "\n")

md5     <- package_source_md5()
reasons <- validate_cached_fit(
  name        = "bcg_meta-analysis",
  deep        = FALSE,
  package_md5 = md5
)
cat("validation reasons for bcg_meta-analysis:\n")
print(reasons)

cat("\nR/zplot-api.R is cache-relevant: ",
    "R/zplot-api.R" %in% .fit_cache_required_source_files(), "\n")
cat("R/fit.R is cache-relevant: ",
    "R/fit.R" %in% .fit_cache_required_source_files(), "\n")
cat("R/utilities.R is cache-relevant: ",
    "R/utilities.R" %in% .fit_cache_required_source_files(), "\n")
