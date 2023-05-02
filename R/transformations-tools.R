#### internal transformation functions ####
### main
# jac = d inv(x) / dx

.d2r <- list(
  fun = function(x) .d2r_fun(x),
  inv = function(x) .r2d$fun(x),
  jac = function(x) 2 / (1 - x^2)^(3/2)
)

.r2d <- list(
  fun = function(x) 2 * x/sqrt(1 - x^2),
  inv = function(x) .d2r$fun(x),
  jac = function(x) 4 / (x^2 + 4)^(3/2)
)

.d2logOR <- list(
  fun = function(x) x * pi/sqrt(3),
  inv = function(x) .logOR2d$fun(x),
  jac = function(x) sqrt(3)/pi
)

.logOR2d <- list(
  fun = function(x) x * sqrt(3)/pi,
  inv = function(x) .d2logOR$fun(x),
  jac = function(x) pi/sqrt(3)
)

.z2r <- list(
  fun = function(x) (exp(2 * x) - 1)/(1 + exp(2 * x)),
  inv = function(x) .r2z$fun(x),
  jac = function(x) - 1 / (x^2 - 1)
)

.r2z <- list(
  fun = function(x) 0.5 * log((1 + x)/(1 - x)),
  inv = function(x) .z2r$fun(x),
  jac = function(x) 4 * exp(2*x) / (exp(2*x) + 1)^2
)

.OR2logOR <- list(
  fun = function(x) log(x),
  inv = function(x) .logOR2OR$fun(x),
  jac = function(x) exp(x)
)

.logOR2OR <- list(
  fun = function(x) exp(x),
  inv = function(x) .OR2logOR$fun(x),
  jac = function(x) 1/x
)

# d2r needs to check for infinity
.d2r_fun <- function(d){
  r <- d/sqrt(d^2 + 4)

  # and possibility of passing data.frames
  if(is.data.frame(r)){
    for(i in 1:ncol(r)){
      r[is.infinite(d[,i]) & d[,i] > 0, i] <-  1
      r[is.infinite(d[,i]) & d[,i] < 0, i] <- -1
    }
  }else{
    r[is.infinite(d) & d > 0] <-  1
    r[is.infinite(d) & d < 0] <- -1
  }

  return(r)
}

### compound
# fun = t2(t1(x))
# jac = d inv(x) / dx = d inv_t2(inv_t1(x)) / dx = jac_t1(inv_t1(x)) * jac_t2 [by chain rule]

.d2z <- list(
  fun = function(x) .r2z$fun(.d2r$fun(x)),
  inv = function(x) .z2d$fun(x),
  jac = function(x) .d2r$jac(.z2r$fun(x)) * .r2z$jac(x)
)

.z2d <- list(
  fun = function(x) .r2d$fun(.z2r$fun(x)),
  inv = function(x) .d2z$fun(x),
  jac = function(x) .z2r$jac(.d2r$fun(x)) * .r2d$jac(x)
)

.logOR2z <- list(
  fun = function(x) .d2z$fun(.logOR2d$fun(x)),
  inv = function(x) .z2logOR$fun(x),
  jac = function(x) .logOR2d$jac(.z2d$fun(x)) * .d2z$jac(x)
)

.z2logOR <- list(
  fun = function(x) .d2logOR$fun(.z2d$fun(x)),
  inv = function(x) .logOR2z$fun(x),
  jac = function(x) .z2d$jac(.logOR2d$fun(x)) * .d2logOR$jac(x)
)

.logOR2r <- list(
  fun = function(x) .d2r$fun(.logOR2d$fun(x)),
  inv = function(x) .r2logOR$fun(x),
  jac = function(x) .logOR2d$jac(.r2d$fun(x)) * .d2r$jac(x)
)

.r2logOR <- list(
  fun = function(x) .d2logOR$fun(.r2d$fun(x)),
  inv = function(x) .logOR2r$fun(x),
  jac = function(x) .r2d$jac(.logOR2d$fun(x)) * .d2logOR$jac(x)
)

.d2OR <- list(
  fun = function(x) .logOR2OR$fun(.d2logOR$fun(x)),
  inv = function(x) .OR2d$fun(x),
  jac = function(x) .d2logOR$jac(.OR2logOR$fun(x)) * .logOR2OR$jac(x)
)

.OR2d <- list(
  fun = function(x) .logOR2d$fun(.OR2logOR$fun(x)),
  inv = function(x) .d2OR$fun(x),
  jac = function(x) .OR2logOR$jac(.d2logOR$fun(x)) * .logOR2d$jac(x)
)

.r2OR <- list(
  fun = function(x) .logOR2OR$fun(.r2logOR$fun(x)),
  inv = function(x) .OR2r$fun(x),
  jac = function(x) .r2logOR$jac(.OR2logOR$fun(x)) * .logOR2OR$jac(x)
)

.OR2r <- list(
  fun = function(x) .logOR2r$fun(.OR2logOR$fun(x)),
  inv = function(x) .r2OR$fun(x),
  jac = function(x) .OR2logOR$jac(.r2logOR$fun(x)) * .logOR2r$jac(x)
)

.z2OR <- list(
  fun = function(x) .logOR2OR$fun(.z2logOR$fun(x)),
  inv = function(x) .OR2z$fun(x),
  jac = function(x) .z2logOR$jac(.OR2logOR$fun(x)) * .logOR2OR$jac(x)
)

.OR2z <- list(
  fun = function(x) .logOR2z$fun(.OR2logOR$fun(x)),
  inv = function(x) .z2OR$fun(x),
  jac = function(x) .OR2logOR$jac(.z2logOR$fun(x)) * .logOR2z$jac(x)
)

#### internal scaling functions ####
### main (could've been simplified but this keeps them parallel with the transformations)
# jac = d inv(x) / dx
# the r/OR scaling is only used for PET/PEESE coefficients and exists for consistency (it's not used for heterogeneity)

.scale_d2r <- list(
  fun = function(x) x / 2,
  inv = function(x) .scale_r2d$fun(x),
  jac = function(x) 2
)

.scale_r2d <- list(
  fun = function(x) 2 * x,
  inv = function(x) .scale_d2r$fun(x),
  jac = function(x) 1/2
)

.scale_d2logOR <- list(
  fun = function(x) x * pi/sqrt(3),
  inv = function(x) .scale_logOR2d$fun(x),
  jac = function(x) sqrt(3)/pi
)

.scale_logOR2d <- list(
  fun = function(x) x * sqrt(3)/pi,
  inv = function(x) .scale_d2logOR$fun(x),
  jac = function(x) pi/sqrt(3)
)

.scale_z2r <- list(
  fun = function(x) x,
  inv = function(x) .scale_r2z$fun(x),
  jac = function(x) 1
)

.scale_r2z <- list(
  fun = function(x) x,
  inv = function(x) .scale_z2r$fun(x),
  jac = function(x) 1
)

.scale_logOR2OR <- list(
  fun = function(x) exp(x),
  inv = function(x) .scale_OR2logOR,
  jac = function(x) 1/x
)

.scale_OR2logOR <- list(
  fun = function(x) log(x),
  inv = function(x) .scale_logOR2OR,
  jac = function(x) exp(x)
)

### compound
# fun = t2(t1(x))
# jac = d inv(x) / dx = d inv_t2(inv_t1(x)) / dx = jac_t1(inv_t1(x)) * jac_t2 [by chain rule]

.scale_d2z <- list(
  fun = function(x) .scale_r2z$fun(.scale_d2r$fun(x)),
  inv = function(x) .scale_z2d$fun(x),
  jac = function(x) .scale_d2r$jac(.scale_z2r$fun(x)) * .scale_r2z$jac(x)
)

.scale_z2d <- list(
  fun = function(x) .scale_r2d$fun(.scale_z2r$fun(x)),
  inv = function(x) .scale_d2z$fun(x),
  jac = function(x) .scale_z2r$jac(.scale_d2r$fun(x)) * .scale_r2d$jac(x)
)

.scale_logOR2z <- list(
  fun = function(x) .scale_d2z$fun(.scale_logOR2d$fun(x)),
  inv = function(x) .scale_z2logOR$fun(x),
  jac = function(x) .scale_logOR2d$jac(.scale_z2d$fun(x)) * .scale_d2z$jac(x)
)

.scale_z2logOR <- list(
  fun = function(x) .scale_d2logOR$fun(.scale_z2d$fun(x)),
  inv = function(x) .scale_logOR2z$fun(x),
  jac = function(x) .scale_z2d$jac(.scale_logOR2d$fun(x)) * .scale_d2logOR$jac(x)
)

.scale_logOR2r <- list(
  fun = function(x) .scale_d2r$fun(.scale_logOR2d$fun(x)),
  inv = function(x) .scale_r2logOR$fun(x),
  jac = function(x) .scale_logOR2d$jac(.scale_r2d$fun(x)) * .scale_d2r$jac(x)
)

.scale_r2logOR <- list(
  fun = function(x) .scale_d2logOR$fun(.scale_r2d$fun(x)),
  inv = function(x) .scale_logOR2r$fun(x),
  jac = function(x) .scale_r2d$jac(.scale_logOR2d$fun(x)) * .scale_d2logOR$jac(x)
)

.scale_z2OR <- list(
  fun = function(x) .scale_logOR2OR$fun(.scale_z2logOR$fun(x)),
  inv = function(x) .scale_OR2z$fun(x),
  jac = function(x) .scale_z2logOR$jac(.scale_OR2logOR$fun(x)) * .scale_logOR2OR$jac(x)
)

.scale_OR2z <- list(
  fun = function(x) .scale_logOR2z$fun(.scale_logOR2OR$fun(x)),
  inv = function(x) .scale_z2OR$fun(x),
  jac = function(x) .scale_logOR2OR$jac(.scale_z2logOR$fun(x)) * .scale_logOR2z$jac(x)
)

.scale_d2OR <- list(
  fun = function(x) .scale_logOR2OR$fun(.scale_d2logOR$fun(x)),
  inv = function(x) .scale_OR2d$fun(x),
  jac = function(x) .scale_d2logOR$jac(.scale_OR2logOR$fun(x)) * .scale_logOR2OR$jac(x)
)

.scale_OR2d <- list(
  fun = function(x) .scale_logOR2d$fun(.scale_logOR2OR$fun(x)),
  inv = function(x) .scale_d2OR$fun(x),
  jac = function(x) .scale_logOR2OR$jac(.scale_d2logOR$fun(x)) * .scale_logOR2d$jac(x)
)

.scale_r2OR <- list(
  fun = function(x) .scale_logOR2OR$fun(.scale_r2logOR$fun(x)),
  inv = function(x) .scale_OR2r$fun(x),
  jac = function(x) .scale_r2logOR$jac(.scale_OR2logOR$fun(x)) * .scale_logOR2OR$jac(x)
)

.scale_OR2r <- list(
  fun = function(x) .scale_logOR2r$fun(.scale_logOR2OR$fun(x)),
  inv = function(x) .scale_r2OR$fun(x),
  jac = function(x) .scale_logOR2OR$jac(.scale_r2logOR$fun(x)) * .scale_logOR2r$jac(x)
)
