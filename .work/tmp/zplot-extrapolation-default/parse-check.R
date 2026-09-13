files <- c(
  "R/zplot-api.R",
  "tests/testthat/test-00-zplot-defaults.R",
  "tests/testthat/test-00-selection-kernel-plots.R",
  "tests/testthat/test-03-zplot.R",
  "tests/scenarios/test-assink2016.R",
  "tests/scenarios/test-hoogeveen2023.R",
  "tests/scenarios/test-ishak2007.R",
  "tests/scenarios/test-kearon1998.R",
  "tests/scenarios/test-white2020.R"
)

for (file in files) {
  parsed <- tryCatch(parse(file), error = function(e) e)
  if (inherits(parsed, "error")) {
    cat("PARSE FAILED:", file, "-", conditionMessage(parsed), "\n")
  } else {
    cat("ok:", file, "(", length(parsed), "expressions )\n")
  }
}
