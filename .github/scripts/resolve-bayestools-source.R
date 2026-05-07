description <- read.dcf("DESCRIPTION")[1, ]
imports     <- description[["Imports"]]

match <- regexec(
  "BayesTools\\s*\\(\\s*>=\\s*([^)[:space:]]+)\\s*\\)",
  imports,
  perl = TRUE
)

required_version <- regmatches(imports, match)[[1]]
required_version <- if (length(required_version) >= 2) required_version[[2]] else NA_character_

cran_version <- tryCatch(
  {
    repos <- getOption("repos")
    if (is.null(repos) ||
        is.na(repos[["CRAN"]]) ||
        identical(unname(repos[["CRAN"]]), "@CRAN@")) {
      repos <- c(CRAN = "https://cloud.r-project.org")
    }

    available <- utils::available.packages(
      contriburl = utils::contrib.url(repos, type = "source")
    )
    available["BayesTools", "Version"]
  },
  error = function(e) NA_character_
)

use_github   <- !is.na(required_version) &&
  (is.na(cran_version) ||
     utils::compareVersion(required_version, cran_version) > 0)

pak_spec     <- if (use_github) {
  "github::FBartos/BayesTools@major-refactoring"
} else {
  "any::BayesTools"
}

message(
  "BayesTools required: ", ifelse(is.na(required_version), "none", required_version),
  "; CRAN: ", ifelse(is.na(cran_version), "unavailable", cran_version),
  "; pak spec: ", pak_spec
)

output_file <- Sys.getenv("GITHUB_OUTPUT")
if (nzchar(output_file)) {
  cat("pak_spec=", pak_spec, "\n", sep = "", file = output_file, append = TRUE)
} else {
  cat(pak_spec, "\n", sep = "")
}
