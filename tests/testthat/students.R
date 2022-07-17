students <- data.frame(readxl::read_excel("C:/Users/fbart/Downloads/students.xlsx"))

head(students)

# create a shorter data-set
df <- data.frame(
  # effect size and standard errors
  r  = students$pcc,
  se = students$se_pcc,

  study       = students$study,
  sample_size = students$sample_size_full,

  # outcome type
  education_outcome = ifelse(students$education_test == 1, "test",
                     ifelse(students$education_choice == 1, "choice",
                            ifelse(students$education_attainment == 1, "attainment", NA))),

  employment_intensity = ifelse(students$low_intensity_employment == 1, "low",
                                ifelse(students$medium_intensity_employment == 1, "medium",
                                       ifelse(students$high_intensity_employment == 1, "high", NA))),

  # students
  students_gender = ifelse(students$male_students == 1, "male",
               ifelse(students$female_students == 1, "female",
                      ifelse(students$mixed_gender_students == 1, "mixed", NA))),

  # location
  location = ifelse(students$usa == 1, "USA",
                    ifelse(students$germany == 1, "Germany",
                           ifelse(students$europe == 1, "Europe", "Other"))),

  # design
  design = ifelse(students$longitudinal_data == 1, "longitudinal",
                  ifelse(students$crosssectional_data == 1, "crosssectional", NA)),

  # endogeneity
  endogenity_control = ifelse(students$endogeneity_control == 1, "Yes", "No"),

  # motivation
  motivation_control = ifelse(students$motivation_control == 1, "Yes", "No")
)

Kroupova2021 <- df
save(Kroupova2021, file = "data/Kroupova2021.RData", version = 2)


length(names(table(df$study)))

df      <- na.omit(df)
df$z_se <- df$se

# reproduce Table 3 (without clustered coefficients)
summary(lm(r ~ 1 + se, data = df, weights = 1/se^2))   # ~ OLS
summary(lm(r ~ 1 + I(1/sqrt(sample_size)), data = df)) # ~ IV

summary(lm(r ~ 1 + I(1/sqrt(sample_size)), data = df))


fit_linear0 <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = df,
  study_ids = df$study,

  # do not test for the presence of the predictors & and do not scale the predictors
  test_predictors        = "",
  standardize_predictors = FALSE,

  # specify very vague priors for the effect and regression coefficient (on Fisher's z scale)
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 10)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_effect_null   = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)

fit_linear <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ z_se,
  data      = df,
  study_ids = df$study,

  # do not test for the presence of the predictors & and do not scale the predictors
  test_predictors        = "",
  standardize_predictors = FALSE,

  # specify very vague priors for the effect and regression coefficient (on Fisher's z scale)
  priors        = list("z_se" = prior("normal", parameters = list(mean = 0, sd = 10))),
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 10)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_effect_null   = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)

fit_linear2 <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ z_se,
  data      = df,

  # do not test for the presence of the predictors & and do not scale the predictors
  test_predictors        = "",
  standardize_predictors = FALSE,

  # specify very vague priors for the effect and regression coefficient (on Fisher's z scale)
  priors        = list("z_se" = prior("normal", parameters = list(mean = 0, sd = 10))),
  priors_effect = prior("normal", parameters = list(mean = 0, sd = 10)),
  prior_scale   = "fishers_z",

  # remove the remaining model components
  priors_bias          = NULL,
  priors_effect_null   = NULL,
  priors_heterogeneity = NULL,

  # some additional settings
  parallel = TRUE, seed = 1
)


fit_linear_random <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ z_se,
  data      = df,
  study_ids = df$study,

  # do not test for the presence of the predictors & and do not scale the predictors
  test_predictors        = "",
  standardize_predictors = FALSE,

  # specify very vague priors for the effect and regression coefficient (on Fisher's z scale)
  priors               = list("z_se" = prior("normal", parameters = list(mean = 0, sd = 10))),
  priors_effect        = prior("normal", parameters = list(mean = 0, sd = 10)),
  priors_heterogeneity = prior("normal", parameters = list(mean = 0, sd = 10), truncation = list(lower = 0)),
  prior_scale          = "fishers_z",

  # remove the remaining model components
  priors_bias               = NULL,
  priors_effect_null        = NULL,
  priors_heterogeneity_null = NULL,

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)

summary(fit_RoBMA, output_scale = "r")
summary(fit_RoBMA, output_scale = "fishers_z")
fit_RoBMA <- RoBMA.reg(
  # specify the model formula and data input
  formula   = ~ 1,
  data      = df,
  study_ids = df$study,

  # specify the default RoBMA model
  model_type       = "PSMA",
  effect_direction = "negative",

  # some additional settings
  parallel = TRUE, seed = 1, weighted = TRUE
)
saveRDS(fit_RoBMA, file = "models/Students/fit_RoBMA.RDS")
saveRDS(fit_RoBMA, file = "models/Students/fit_RoBMA.RDS")

summary(fit_linear$models[[1]]$fit)
summary(fit_linear2$models[[1]]$fit)

summary(fit_linear, type = "individual", output_scale = "r")
summary(fit_linear2, type = "individual", output_scale = "r")
summary(fit_linear_random, type = "individual", output_scale = "r")

summary(fit_linear, type = "individual", output_scale = "fishers_z")
summary(fit_linear_random, type = "individual", output_scale = "fishers_z")

summary(fit_fe1w,output_scale = "fishers_z", )

fit_metafor <- metafor::rma(
  yi     = df$z,
  sei    = df$se^2,
  mods   = ~ education + location + motivation,
  data   = df,
)
fit_metafor <- metafor::rma.mv(
  yi     = df$z,
  V      = df$se^2,
  mods   = ~ educfation + location + motivation,
  random = ~ 1 | study,
  data   = df
)
summary(fit_metafor)
saveRDS(fit_metafor, file = "fit_metafor.RDS")

object$models [[1]]$fit_summary
object$models [[1]]$fit_summaries$z

fit_normal <- RoBMA.reg(
  formula          = ~ education + location + motivation,
  test_predictors  = c("education", "location", "motivation"),
  data             = df,
  study_ids        = df$study,
  effect_direction = "negative",

  # no publication bias adjustment
  priors_bias     = NULL,

  # fairly tight priors for the BF-tests of the covariates
  prior_factors   = prior_factor("mnorm", parameters = list(mean = 0, sd = 0.25)),

  parallel        = TRUE,
  weighted        = TRUE
)
saveRDS(fit_normal, file = "fit_normal.RDS")

fit_adjusted <- RoBMA.reg(
  formula          = ~ education + location + motivation,
  test_predictors  = c("education", "location", "motivation"),
  data             = df,
  study_ids        = df$study,
  effect_direction = "negative",

  # fairly tight priors for the BF-tests of the covariates
  prior_factors   = prior_factor("mnorm", parameters = list(mean = 0, sd = 0.25)),

  parallel        = TRUE,
  chains          = 1,
  sample          = 500,
  burnin          = 200,
  weighted        = TRUE
)
saveRDS(fit_adjusted, file = "fit_adjusted.RDS")

z2r(students$estimate)
students$estimate
students$se_estimate


datagrades = read.table("clipboard-512", sep="\t", header=TRUE)
grades0 = bms(datagrades, burn=1e5,iter=3e5, g="UIP", mprior="uniform", nmodel=50000, mcmc="bd", user.int=FALSE)
grades = bms(datagrades, burn=1e5,iter=3e5, g="UIP", mprior="dilut", nmodel=50000, mcmc="bd", user.int=FALSE)
grades1 = bms(datagrades, burn=1e5,iter=3e5, g="BRIC", mprior="random", nmodel=50000, mcmc="bd", user.int=FALSE)
coef(grades, order.by.pip = F, exact=T, include.constant=T)
image(grades, yprop2pip=FALSE, order.by.pip=TRUE, do.par=TRUE, do.grid=TRUE, do.axis=TRUE, cex.axis = 0.7)
summary(grades)
plot(grades)
print(grades$topmod[1])


par(mfrow=c(5,2))
density(grades, reg="se_no_endogeneity_control")
density(grades, reg="ols_method")
density(grades, reg="germany")
density(grades, reg="longitudinal_data")
density(grades, reg="education_choice")
density(grades, reg="employment_continuous")
density(grades, reg="low_intensity_employment")
density(grades, reg="high_intensity_employment")
density(grades, reg="ability_control")
density(grades, reg="ethnicity_control")

*drop <- c("data_year")

library(corrplot)
datagrades = read.table("clipboard-512", sep="\t", header=TRUE)
col<- colorRampPalette(c("red", "white", "blue"))
M <- cor(datagrades)
corrplot.mixed(M, lower = "number", upper = "circle", lower.col=col(200), upper.col=col(200), tl.pos = c("lt"), diag = c("u"), tl.col="black", tl.srt=45, tl.cex=0.85, number.cex = 0.5, cl.cex=0.8, cl.ratio=0.1)


library(foreign)
library(xtable)
library(LowRankQP)
datagrades=read.table("clipboard-512", sep="\t", header=TRUE)
datagrades <-na.omit(datagrades)
x.data <- datagrades[,-1]
const_<-c(1)
x.data <-cbind(const_,x.data)

x <- sapply(1:ncol(x.data),function(i){x.data[,i]/max(x.data[,i])})
scale.vector <- as.matrix(sapply(1:ncol(x.data),function(i){max(x.data[,i])}))
Y <- as.matrix(datagrades[,1])
output.colnames <- colnames(x.data)
full.fit <- lm(Y~x-1)
beta.full <- as.matrix(coef(full.fit))
M <- k <- ncol(x)
n <- nrow(x)
beta <- matrix(0,k,M)
e <- matrix(0,n,M)
K_vector <- matrix(c(1:M))
var.matrix <- matrix(0,k,M)
bias.sq <- matrix(0,k,M)

for(i in 1:M)
{
  X <- as.matrix(x[,1:i])
  ortho <- eigen(t(X)%*%X)
  Q <- ortho$vectors ; lambda <- ortho$values
  x.tilda <- X%*%Q%*%(diag(lambda^-0.5,i,i))
  beta.star <- t(x.tilda)%*%Y
  beta.hat <- Q%*%diag(lambda^-0.5,i,i)%*%beta.star
  beta[1:i,i] <- beta.hat
  e[,i] <- Y-x.tilda%*%as.matrix(beta.star)
  bias.sq[,i] <- (beta[,i]-beta.full)^2
  var.matrix.star <- diag(as.numeric(((t(e[,i])%*%e[,i])/(n-i))),i,i)
  var.matrix.hat <- var.matrix.star%*%(Q%*%diag(lambda^-1,i,i)%*%t(Q))
  var.matrix[1:i,i] <- diag(var.matrix.hat)
  var.matrix[,i] <- var.matrix[,i]+ bias.sq[,i]
}

e_k <- e[,M]
sigma_hat <- as.numeric((t(e_k)%*%e_k)/(n-M))
G <- t(e)%*%e
a <- ((sigma_hat)^2)*K_vector
A <- matrix(1,1,M)
b <- matrix(1,1,1)
u <- matrix(1,M,1)
optim <- LowRankQP(Vmat=G,dvec=a,Amat=A,bvec=b,uvec=u,method="LU",verbose=FALSE)
weights <- as.matrix(optim$alpha)
beta.scaled <- beta%*%weights
final.beta <- beta.scaled/scale.vector
std.scaled <- sqrt(var.matrix)%*%weights
final.std <- std.scaled/scale.vector
results.reduced <- as.matrix(cbind(final.beta,final.std))
rownames(results.reduced) <- output.colnames; colnames(results.reduced) <- c("Coefficient", "Sd. Err")
MMA.fls <- round(results.reduced,4)
MMA.fls <- data.frame(MMA.fls)
t <- as.data.frame(MMA.fls$Coefficient/MMA.fls$Sd..Err)
MMA.fls$pv <-round( (1-apply(as.data.frame(apply(t,1,abs)), 1, pnorm))*2,3)
MMA.fls$names <- rownames(MMA.fls)
names <- c(colnames(datagrades))
names <- c(names,"const_")
MMA.fls <- MMA.fls[match(names, MMA.fls$names),]
MMA.fls$names <- NULL
MMA.fls
