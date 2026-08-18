#' @title 23 experimental studies from
#' \insertCite{anderson2010violent;textual}{RoBMA} that meet the best practice criteria
#'
#' @description The data set contains correlation coefficients, sample
#' sizes, and labels for 23 experimental studies focusing on the effect of
#' violent video games on aggressive behavior.
#'
#'
#' @format A data.frame with 3 columns and 23 observations:
#' \describe{
#'   \item{\code{r}}{Correlation coefficient.}
#'   \item{\code{n}}{Sample size.}
#'   \item{\code{name}}{Study label.}
#' }
#' @source \url{https://github.com/Joe-Hilgard/Anderson-meta}
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
#' @format A data.frame with 3 columns and 9 observations:
#' \describe{
#'   \item{\code{d}}{Cohen's d effect size.}
#'   \item{\code{se}}{Standard error of d.}
#'   \item{\code{study}}{Study label.}
#' }
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
#' @format A data.frame with 3 columns and 5 observations:
#' \describe{
#'   \item{\code{d}}{Cohen's d effect size.}
#'   \item{\code{se}}{Standard error of d.}
#'   \item{\code{study}}{Study label.}
#' }
#' @references
#' \insertAllCited{}
"Poulsen2006"

#' @title 881 estimates from 69 studies of a relationship between employment and
#' educational outcomes collected by \insertCite{kroupova2021student;textual}{RoBMA}
#'
#' @description The data set contains partial correlation coefficients, standard errors,
#' study labels, sample sizes, type of the educational outcome, intensity of the
#' employment, gender of the student population, study location, study design, whether
#' the study controlled for endogeneity, and whether the study controlled for motivation.
#' The original data set including additional variables and the publication are
#' available from the project page.
#' (Note that some standard errors and employment intensities are missing.)
#'
#' @format A data.frame with 11 columns and 881 observations:
#' \describe{
#'   \item{\code{r}}{Partial correlation coefficient.}
#'   \item{\code{se}}{Standard error of r.}
#'   \item{\code{study}}{Study label.}
#'   \item{\code{sample_size}}{Sample size.}
#'   \item{\code{education_outcome}}{Type of educational outcome.}
#'   \item{\code{employment_intensity}}{Employment intensity.}
#'   \item{\code{students_gender}}{Gender composition of the student sample.}
#'   \item{\code{location}}{Study location.}
#'   \item{\code{design}}{Study design.}
#'   \item{\code{endogenity_control}}{Whether endogeneity was controlled.}
#'   \item{\code{motivation_control}}{Whether motivation was controlled.}
#' }
#' @source \url{https://meta-analysis.cz/students/}
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
#' @format A data.frame with 3 columns and 18 observations:
#' \describe{
#'   \item{\code{r}}{Correlation coefficient.}
#'   \item{\code{n}}{Sample size.}
#'   \item{\code{study}}{Study label.}
#' }
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
#' @format A data.frame with 12 columns and 70 observations:
#' \describe{
#'   \item{\code{Learning_effect}}{Learning-effect outcome category.}
#'   \item{\code{Author_year}}{Author-year study label.}
#'   \item{\code{N_EG}}{Experimental-group sample size.}
#'   \item{\code{N_CG}}{Control-group sample size.}
#'   \item{\code{g}}{Hedges' g effect size.}
#'   \item{\code{Grade_level}}{Grade level.}
#'   \item{\code{Type_of_course}}{Course type.}
#'   \item{\code{Duration}}{Study duration.}
#'   \item{\code{Learning_model}}{Learning model.}
#'   \item{\code{Role_of_ChatGPT}}{Role of ChatGPT in the intervention.}
#'   \item{\code{Area_of_ChatGPT_application}}{Application area.}
#'   \item{\code{se}}{Standard error of g.}
#' }
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
#' @format A data.frame with 2 columns and 55 observations:
#' \describe{
#'   \item{\code{y}}{Standardized mean difference.}
#'   \item{\code{se}}{Standard error of y.}
#' }
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
#' @format A data.frame with 9 columns and 37 observations:
#' \describe{
#'   \item{\code{d}}{Cohen's d effect size.}
#'   \item{\code{v}}{Sampling variance of d.}
#'   \item{\code{outcome}}{Outcome type.}
#'   \item{\code{feedback_level}}{Feedback level.}
#'   \item{\code{social_comparison_type}}{Social-comparison type.}
#'   \item{\code{sessions}}{Number of sessions.}
#'   \item{\code{sample_type}}{Sample type.}
#'   \item{\code{sample_size}}{Sample size.}
#'   \item{\code{country}}{Country.}
#' }
#' @references
#' \insertAllCited{}
"Hoppen2025"

#' @title 582 effect sizes examining the ease-of-retrieval effect
#' from a meta-analysis by \insertCite{weingarten2018does;textual}{RoBMA}
#'
#' @description The data set contains correlation coefficients between the manipulation
#' and outcome variable (r_xy), sample sizes (N), and various study characteristics
#' including publication status, country (USA vs. other), number of few and many
#' examples requested, whether the trial targeted episodic memory, paradigm type (standard vs. other),
#' dataset type (proximal vs. distal), and mediation variables (r_xm, r_my).
#' The meta-analysis examined whether subjective ease mediates the ease-of-retrieval effect,
#' where participants list either few or many examples and then make judgments
#' \insertCite{weingarten2018does}{RoBMA}.
#'
#' @format A data.frame with 12 columns and 582 observations:
#' \describe{
#'   \item{\code{r_xy}}{Correlation between manipulation and outcome.}
#'   \item{\code{N}}{Sample size.}
#'   \item{\code{paper_id}}{Paper identifier.}
#'   \item{\code{published}}{Publication-status indicator.}
#'   \item{\code{USA}}{Whether the study was conducted in the USA.}
#'   \item{\code{number_of_few}}{Number of examples in the few-examples condition.}
#'   \item{\code{number_of_many}}{Number of examples in the many-examples condition.}
#'   \item{\code{episodic_memory}}{Whether the trial targeted episodic memory.}
#'   \item{\code{standard_paradigm}}{Whether the standard paradigm was used.}
#'   \item{\code{proximal_dataset}}{Whether the data set was proximal.}
#'   \item{\code{r_xm}}{Correlation between manipulation and mediator.}
#'   \item{\code{r_my}}{Correlation between mediator and outcome.}
#' }
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
#' @format A data.frame with 3 columns and 412 observations:
#' \describe{
#'   \item{\code{study}}{Study label.}
#'   \item{\code{d}}{Cohen's d effect size.}
#'   \item{\code{se}}{Standard error of d.}
#' }
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
#' @format A data.frame with 5 columns and 1159 observations:
#' \describe{
#'   \item{\code{y}}{Effect size: percent increase in earnings.}
#'   \item{\code{se}}{Standard error of y.}
#'   \item{\code{facing_customer}}{Type of customer contact.}
#'   \item{\code{study_id}}{Study identifier.}
#'   \item{\code{N}}{Sample size.}
#' }
#' @references
#' \insertAllCited{}
"Havrankova2025"


#' @title 106 analyst-team effect-size estimates of the relation between
#' religiosity and well-being from \insertCite{hoogeveen2023many;textual}{RoBMA}
#'
#' @description The data set contains standardized effect-size estimates,
#' standard errors, 95\% confidence interval bounds, sample sizes, effect-size
#' types, sequential analyst-team identifiers, and team-level characteristics
#' describing psychological dependent-variable coding and self-reported
#' methodological knowledge. The Many-Analysts Religion Project examined
#' whether religious people self-report greater well-being
#' \insertCite{hoogeveen2023many}{RoBMA}.
#'
#' @format A data.frame with 9 columns and 106 observations:
#' \describe{
#'   \item{\code{yi}}{Standardized effect-size estimate.}
#'   \item{\code{sei}}{Standard error of the effect-size estimate.}
#'   \item{\code{lci}}{Lower bound of the reported 95\% confidence interval.}
#'   \item{\code{uci}}{Upper bound of the reported 95\% confidence interval.}
#'   \item{\code{ni}}{Analysis sample size.}
#'   \item{\code{type}}{Reported effect-size type: standardized regression coefficient beta or Cohen's d.}
#'   \item{\code{team}}{Sequential analyst-team identifier in this example.}
#'   \item{\code{team_psych}}{Whether the effect estimate used the psychological dependent-variable coding.}
#'   \item{\code{team_knowledge}}{Team's self-reported methodological-knowledge rating, from 1 (no knowledge) to 5 (expert).}
#' }
#' @source \url{https://osf.io/gkxqy/}
#'
#' @references
#' \insertAllCited{}
"Hoogeveen2023"
