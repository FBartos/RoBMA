library(RoBMA)
x1 <- c(0, 5, 10,  5, 3)
x2 <- c(1, 4, 15, 10, 8)
n1 <- c(8, 15,30, 30,10)
n2 <- c(7, 14,32, 30,10)

debug(RoBMA:::.fit_BiBMA_model)
fit <- BiBMA(x1 = x1, x2 = x2, n1 = n1, n2 = n2, silent = FALSE, seed = 1)
