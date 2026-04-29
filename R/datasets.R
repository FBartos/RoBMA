#' @title 27 experimental studies from
#' \insertCite{anderson2010violent;textual}{RoBMA} that meet the best practice criteria
#'
#' @description The data set contains correlation coefficients, sample
#' sizes, and labels for 27 experimental studies focusing on the effect of
#' violent video games on aggressive behavior. The full original data can
#' found at https://github.com/Joe-Hilgard/Anderson-meta.
#'
#'
#' @format A data.frame with 3 columns and 23 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Anderson2010"


#' @title 9 experimental studies from
#' \insertCite{bem2011feeling;textual}{RoBMA} as described in
#' \insertCite{bem2011must;textual}{RoBMA}
#'
#' @description The data set contains Cohen's d effect sizes, standard errors,
#' and labels for 9 experimental studies of precognition from the infamous
#' \insertCite{bem2011feeling;textual}{RoBMA} as analyzed in his later meta-analysis
#' \insertCite{bem2011must}{RoBMA}.
#'
#' @format A data.frame with 3 columns and 9 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Bem2011"


#' @title 5 studies with a tactile outcome assessment from
#' \insertCite{poulsen2006potassium;textual}{RoBMA} of the effect of potassium-containing toothpaste
#' on dentine hypersensitivity
#'
#' @description The data set contains Cohen's d effect sizes, standard errors,
#' and labels for 5 studies assessing the tactile outcome from a meta-analysis of
#' the effect of potassium-containing toothpaste on dentine hypersensitivity
#' \insertCite{poulsen2006potassium}{RoBMA} which was used as an example in
#' \insertCite{bartos2021bayesian;textual}{RoBMA}.
#'
#' @format A data.frame with 3 columns and 5 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Poulsen2006"

#' @title 881 estimates from 69 studies of a relationship between employment and
#' educational outcomes collected by \insertCite{kroupova2021student;textual}{RoBMA}
#'
#' @description The data set contains partial correlation coefficients, standard errors,
#' study labels, samples sizes, type of the educational outcome, intensity of the
#' employment, gender of the student population, study location, study design, whether
#' the study controlled for endogenity, and whether the study controlled for motivation.
#' The original data set including additional variables and the publication can be found
#' at http://meta-analysis.cz/students.
#' (Note that some standard errors and employment intensities are missing.)
#'
#' @format A data.frame with 11 columns and 881 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Kroupova2021"

#' @title 18 studies of a relationship between acculturation mismatch and
#' intergenerational cultural conflict collected by
#' \insertCite{lui2015intergenerational;textual}{RoBMA}
#'
#' @description The data set contains correlation coefficients r,
#' sample sizes n, and labels for each study assessing the
#' relationship between acculturation mismatch (that is the result of the contrast
#' between the collectivist cultures of Asian and Latin immigrant groups
#' and the individualist culture in the United States) and intergenerational cultural
#' conflict \insertCite{lui2015intergenerational}{RoBMA} which was used as an
#' example in \insertCite{bartos2020adjusting;textual}{RoBMA}.
#'
#' @format A data.frame with 3 columns and 18 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Lui2015"

#' @title 39 study rows on household chaos and child executive functions
#' from a meta-analysis by
#' \insertCite{andrews2021examining;textual}{RoBMA}
#'
#' @description The data set contains raw correlation coefficients and
#' sampling information using \pkg{metafor}-compatible column names
#' (\code{ri}, \code{vi}, \code{sei}, \code{ni}, and \code{slab}), together with
#' study-level moderators from the original data set. Three source rows have
#' missing effect-size inputs and are omitted automatically by model-fitting
#' functions that require complete outcomes. The original data set assessed the
#' effect of household chaos on child executive functions
#' \insertCite{andrews2021examining}{RoBMA} and was used as an example in
#' \insertCite{bartos2023robust;textual}{RoBMA}.
#'
#' @format A data.frame with 15 columns and 39 observations:
#' \describe{
#'   \item{\code{study_id}}{Source row identifier.}
#'   \item{\code{slab}}{Study label.}
#'   \item{\code{ri}}{Raw correlation coefficient.}
#'   \item{\code{vi}}{Sampling variance of the raw correlation coefficient.}
#'   \item{\code{sei}}{Sampling standard error of the raw correlation coefficient.}
#'   \item{\code{ni}}{Sample size.}
#'   \item{\code{year}}{Publication year.}
#'   \item{\code{percent_female}}{Percentage of female children in the sample.}
#'   \item{\code{age}}{Mean age of the children in years.}
#'   \item{\code{assessment_interval_months}}{Time between household-chaos and executive-function assessments in months.}
#'   \item{\code{dissertation}}{Whether the study was a dissertation.}
#'   \item{\code{percent_minority}}{Percentage of minority participants in the sample.}
#'   \item{\code{percent_low_parental_education}}{Percentage of parents with high school, GED, or lower education.}
#'   \item{\code{hc_dimension}}{Household-chaos dimension using the source coding.}
#'   \item{\code{measure}}{Executive-function assessment type, direct or informant.}
#' }
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Andrews2021"

#' @title 70 effect sizes from a meta-analysis of ChatGPT's impact on student learning
#' by \insertCite{wang2025effect;textual}{RoBMA}
#'
#' @description The data set contains Hedges' g effect sizes, standard errors,
#' sample sizes for experimental and control groups, and various study characteristics
#' including grade level, type of course, duration, learning model, role of ChatGPT,
#' and area of ChatGPT application. The meta-analysis examined the effect of ChatGPT
#' on students' learning performance, learning perception, and higher-order thinking
#' \insertCite{wang2025effect}{RoBMA}.
#'
#' @format A data.frame with 12 columns and 70 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Wang2025"

#' @title 55 effect sizes from Many Labs 2 replication studies of
#' \insertCite{tversky1981framing;textual}{RoBMA} framing effects
#'
#' @description The data set contains standardized mean differences (y) and standard errors (se)
#' from 55 replication studies of the framing effect originally described by
#' \insertCite{tversky1981framing;textual}{RoBMA}. These studies were part of the
#' Many Labs 2 project examining variation in replicability across samples and settings
#' \insertCite{klein2018many;textual}{RoBMA}.
#'
#' @format A data.frame with 2 columns and 55 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"ManyLabs16"

#' @title 37 studies from a meta-analysis of social comparison as a behavior change technique
#' by \insertCite{hoppen2025meta;textual}{RoBMA}
#'
#' @description The data set contains Cohen's d effect sizes, variances, and study
#' characteristics including outcome type, feedback level, social comparison type,
#' number of sessions, sample type, sample size, and country. The meta-analysis
#' examined social comparison as a behavior change technique across the behavioral
#' sciences \insertCite{hoppen2025meta}{RoBMA}.
#'
#' @format A data.frame with 9 columns and 37 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Hoppen2025"

#' @title 582 effect sizes examining the ease-of-retrieval effect
#' from a meta-analysis by \insertCite{weingarten2018does;textual}{RoBMA}
#'
#' @description The data set contains correlation coefficients between the manipulation
#' and outcome variable (r_xy), sample sizes (N), and various study characteristics
#' including publication status, country (USA vs. other), number of few and many
#' examples requested, whether memory the trial targeted episodic memory, paradigm type (standard vs. other),
#' dataset type (proximal vs. distal), and mediation variables (r_xm, r_my).
#' The meta-analysis examined whether subjective ease mediates the ease-of-retrieval effect,
#' where participants list either few or many examples and then make judgments
#' \insertCite{weingarten2018does}{RoBMA}.
#'
#' @format A data.frame with 12 columns and 582 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Weingarten2018"


#' @title 412 effect sizes from a meta-analysis of secondary benefits of family-based treatments
#' by \insertCite{johnides2025secondary;textual}{RoBMA}
#'
#' @description The data set contains Cohen's d effect sizes, standard errors, and study labels
#' from a meta-analysis investigating the extent to which family-based treatments for children
#' with mental health, physical health, and developmental disorders provide secondary benefits
#' to the children's siblings and caregivers \insertCite{johnides2025secondary}{RoBMA}.
#'
#' @format A data.frame with 3 columns and 412 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Johnides2025"


#' @title 1159 effect sizes from a meta-analysis of beauty and professional success
#' by \insertCite{havrankova2025beauty;textual}{RoBMA}
#'
#' @description The data set contains effect sizes (percent increase in earnings), standard errors,
#' study identifiers, sample sizes, and the type of customer contact (no, some, or direct).
#' The meta-analysis examined the relationship between perceived beauty and professional success
#' \insertCite{havrankova2025beauty}{RoBMA}.
#'
#' @format A data.frame with 5 columns and 1159 observations.
#'
#' @return a data.frame.
#'
#' @references
#' \insertAllCited{}
"Havrankova2025"
