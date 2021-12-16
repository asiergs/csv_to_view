rm(list = ls())

# 1. Data import -------------------------------------------------------------

frec <- read.csv2("https://raw.githubusercontent.com/asiergs/csv_to_view/main/Frecuencia_PF.csv")
frec <- as.numeric(frec[,1])

sev <- read.csv2("https://raw.githubusercontent.com/asiergs/csv_to_view/main/Severidad_PF.csv")
sev <- as.numeric(sev[,1])

mort <- read.csv2("https://raw.githubusercontent.com/asiergs/csv_to_view/main/Mortalidad_PF.csv")

# 2. Estimation method validation (DGP) --------------------------------------

## 2.1  Poisson --------------------------------------------------------------

### 2.1.1 Error functions ----------------------------------------------------

MM_pois <- function(param, sample){
  (mean(sample)-param)^2
}

PM_pois <- function(param, sample){
  (quantile(sample, 0.8)-qpois(0.8, param))^2
}

ML_pois <- function(param, sample){
  -sum(dpois(sample, param, log = TRUE))
}

### 2.1.2 Parameters calculation ---------------------------------------------

n <- length(frec)

param <- 3
sample_DGP <- rpois(n, param)
mean(sample_DGP)

optim(c(1), MM_pois, method = "L-BFGS-B", sample = sample_DGP)
optim(c(1), PM_pois, method = "L-BFGS-B", sample = sample_DGP)
optim(c(1), ML_pois, method = "L-BFGS-B", sample = sample_DGP)

# Since the distribution is discrete, the numerical methods do not properly 
# work. This leads us to the need of another numerical method to estimate the
# parameter based in the Percentile Matching. The following optimization
# function is proposed as a starting point.

optim_discrete <- function(param, FUN,sample, xi=0.1, xf = 20, dx = 1E-1,
                           yi = 0.1, yf = 0.99, dy = 1E-2){
  # for one parameter
  if (length(param)==1){
  range <- seq(xi,xf,dx)
  out <- c()
  for (i in 1:length(range)){
    out[i] <- FUN(range[i],sample)
  }
  sol <- range[out==min(out)]
  list(par = mean(sol), value = min(out))
  }
  # for two parameters
  else if (length(param)==2){
  range_1 <- seq(xi,xf,dx)
  range_2 <- seq(yi,yf,dy)
  out <- matrix(0,length(range_1),length(range_2),
                dimnames = list(range_1, range_2))
  for (i in 1:length(range_1)){
    for (j in 1:length(range_2)){
      param <- c(range_1[i], range_2[j])
      out[i,j] <- FUN(param, sample)
    }
  }
  pos <- which(out == min(out), arr.ind = T)[1,]
  param <- c(range_1[pos[1]], range_2[pos[2]])
  list(par = param, value = min(out))
  }
}
optim_discrete(c(1), PM_pois, sample_DGP)

### 2.1.3 Bias calculation ---------------------------------------------------

MM <- c()
PM <- c()
ML <- c()
sim_size <- 500
n <- length(frec)
param <- 3
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rpois(n, 3)
  MM[i] <- optim(c(1), MM_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  PM[i] <- optim_discrete(c(1), PM_pois, sample_DGP)$par
  ML[i] <- optim(c(1), ML_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  setTxtProgressBar(bar, i)
}
bias_MM <- mean(MM)-param
bias_PM <- mean(PM)-param
bias_ML <- mean(ML)-param

results_pois <- data.frame(method = c("MM", "PM", "ML"),
                              bias = c(bias_MM, bias_PM, bias_ML))

### 2.1.4 MSE calculation ----------------------------------------------------

MSE_MM <- var(MM)+bias_MM^2
MSE_PM <- var(PM)+bias_PM^2
MSE_ML <- var(ML)+bias_ML^2

results_pois <- cbind(results_pois, MSE = c(MSE_MM, MSE_PM, MSE_ML))

### 2.1.5 Consistency --------------------------------------------------------

size_min <- 50
size_max <- 5000
step <- 1
size <- seq(size_min, size_max, step)

MM <- c()
PM <- c()
ML <- c()
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rpois(size[i], param)
  MM[i] <- optim(c(1), MM_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  PM[i] <- optim_discrete(c(1), PM_pois, sample_DGP)$par
  ML[i] <- optim(c(1), ML_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  setTxtProgressBar(bar, i)
}

MM <- MM-param
PM <- PM-param
ML <- ML-param

plot(size, PM, type = "l")
abline(h = 0, col = "red")
# PM method is not consistent
par(mfrow=c(2,1))
plot(size, MM, type = "l")
abline(h = 0, col = "red")
plot(size, ML, type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

## 2.2  Negative Binomial ----------------------------------------------------

### 2.2.1 Error functions ----------------------------------------------------

MM_nbinom <- function(param, sample){
  r1 <- (mean(sample) - param[1]*(1-param[2])/param[2])^2
  r2 <- (var(sample) - param[1]*(1-param[2])/param[2]^2)^2
  r1 + r2
}

PM_nbinom <- function(param, sample){
  r1 <- (quantile(sample, 0.8)-qnbinom(0.8, param[1], param[2]))^2
  r2 <- (quantile(sample, 0.2)-qnbinom(0.2, param[1], param[2]))^2
  r1 + r2
}

ML_nbinom <- function(param, sample){
  -sum(dnbinom(sample, param[1], param[2], log = TRUE))
}

### 2.2.2 Parameters calculation ---------------------------------------------

n <- length(frec)

param <- c(10,0.3)
sample_DGP <- rnbinom(n, param[1], param[2])
mean(sample_DGP)
var(sample_DGP)

optim(c(9,0.4), MM_nbinom, method = "L-BFGS-B", sample = sample_DGP)
optim(c(5,0.5), PM_nbinom, method = "L-BFGS-B", sample = sample_DGP)
optim(c(9,0.4), ML_nbinom, method = "L-BFGS-B", sample = sample_DGP,
      upper = c(NA,0.9999), lower = c(0.00001,0.0001))

optim_discrete(c(1,1), PM_nbinom, sample_DGP)

### 2.2.3 Bias calculation ---------------------------------------------------

sim_size <- 500
MM <- matrix(0,sim_size,2,dimnames = list(c(),c("r", "p")))
PM <- matrix(0,sim_size,2,dimnames = list(c(),c("r", "p")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("r", "p")))

n <- length(frec)
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rnbinom(n, param[1], param[2])
  MM[i,] <- optim(c(9,0.4), MM_nbinom, method = "L-BFGS-B",
                 sample = sample_DGP, lower = c(0.0001,0.0001))$par
  PM[i,] <- optim_discrete(c(5,0.5), PM_nbinom, sample_DGP,
                          xi = 3, xf = 15, dx=0.5, yi = 0.1, yf = 0.5,
                          dy = 0.5E-1)$par
  ML[i,] <- optim(c(9,0.4), ML_nbinom, method = "L-BFGS-B",
                 sample = sample_DGP, upper = c(NA,0.9999),
                 lower = c(0.0001,0.0001))$par
  setTxtProgressBar(bar, i)
}
bias_r_MM <- mean(MM[,1])-param[1]
bias_r_PM <- mean(PM[,1])-param[1]
bias_r_ML <- mean(ML[,1])-param[1]

bias_p_MM <- mean(MM[,2])-param[2]
bias_p_PM <- mean(PM[,2])-param[2]
bias_p_ML <- mean(ML[,2])-param[2]

results_nbinom <- data.frame(method = c("MM", "PM", "ML"),
                             bias_r = c(bias_r_MM, bias_r_PM, bias_r_ML),
                             bias_p = c(bias_p_MM, bias_p_PM, bias_p_ML))

### 2.2.4 MSE calculation ----------------------------------------------------

MSE_r_MM <- var(MM[,1])+bias_r_MM^2
MSE_r_PM <- var(PM[,1])+bias_r_PM^2
MSE_r_ML <- var(ML[,1])+bias_r_ML^2

MSE_p_MM <- var(MM[,2])+bias_p_MM^2
MSE_p_PM <- var(PM[,2])+bias_p_PM^2
MSE_p_ML <- var(ML[,2])+bias_p_ML^2

results_nbinom <- cbind(results_nbinom,
                        MSE_r = c(MSE_r_MM, MSE_r_PM, MSE_r_ML),
                        MSE_p = c(MSE_p_MM, MSE_p_PM, MSE_p_ML))
results_nbinom

### 2.2.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 5000
step <- 3
size <- seq(size_min, size_max, step)

MM <- matrix(0,length(size),2,dimnames = list(c(),c("r", "p")))
PM <- matrix(0,length(size),2,dimnames = list(c(),c("r", "p")))
ML <- matrix(0,length(size),2,dimnames = list(c(),c("r", "p")))
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rnbinom(size[i], param[1], param[2])
  MM[i,] <- optim(c(6,0.6), MM_nbinom, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.0001,0.0001))$par
  PM[i,] <- optim_discrete(c(5,0.5), PM_nbinom, sample_DGP,
                           xi = 3, xf = 15, dx=0.5, yi = 0.1, yf = 0.5,
                           dy = 1E-1)$par
  ML[i,] <- optim(c(6,0.6), ML_nbinom, method = "L-BFGS-B",
                  sample = sample_DGP, upper = c(NA,0.9999),
                  lower = c(0.0001,0.0001))$par
  setTxtProgressBar(bar, i)
}

#### 2.2.5.1 r parameter -----------------------------------------------------

MM[,1] <- MM[,1]-param[1]
PM[,1] <- PM[,1]-param[1]
ML[,1] <- ML[,1]-param[1]

plot(size, PM[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,1], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

#### 2.2.5.2 p parameter -----------------------------------------------------

MM[,2] <- MM[,2]-param[2]
PM[,2] <- PM[,2]-param[2]
ML[,2] <- ML[,2]-param[2]

plot(size, PM[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,2], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

## 2.3 Normal ----------------------------------------------------------------

### 2.3.1 Error functions ----------------------------------------------------

MM_normal<- function(param, sample){
  r1 <- (mean(sample) - param[1])^2
  r2 <- (sd(sample) - param[2])^2
  r1 + r2
}

PM_normal <- function(param, sample){
  r1 <- (quantile(sample, 0.8)-qnorm(0.8, param[1], param[2]))^2
  r2 <- (quantile(sample, 0.2)-qnorm(0.2, param[1], param[2]))^2
  r1 + r2
}

ML_normal <- function(param, sample){
  -sum(dnorm(sample, param[1], param[2], log = TRUE))
}

### 2.3.2 Parameters calculation ---------------------------------------------

n <- length(sev)

param <- c(5000,5000)
sample_DGP <- rnorm(n, param[1], param[2])
mean(sample_DGP)
sd(sample_DGP)

optim(c(4000,4000), MM_normal, method = "L-BFGS-B", sample = sample_DGP)
optim(c(4000,4000), PM_normal, method = "L-BFGS-B", sample = sample_DGP)
optim(c(4000,4000), ML_normal, method = "L-BFGS-B", sample = sample_DGP)

### 2.3.3 Bias calculation ---------------------------------------------------

sim_size <- 500
MM <- matrix(0,sim_size,2,dimnames = list(c(),c("mu", "sigma")))
PM <- matrix(0,sim_size,2,dimnames = list(c(),c("mu", "sigma")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("mu", "sigma")))

n <- length(sev)
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rnorm(n, param[1], param[2])
  MM[i,] <- optim(c(4000,4000), MM_normal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  PM[i,] <- optim(c(4000,4000), PM_normal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  ML[i,] <- optim(c(4000,4000), ML_normal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  setTxtProgressBar(bar, i)
}
bias_mu_MM <- mean(MM[,1])-param[1]
bias_mu_PM <- mean(PM[,1])-param[1]
bias_mu_ML <- mean(ML[,1])-param[1]

bias_sigma_MM <- mean(MM[,2])-param[2]
bias_sigma_PM <- mean(PM[,2])-param[2]
bias_sigma_ML <- mean(ML[,2])-param[2]

results_normal <- data.frame(method = c("MM", "PM", "ML"),
                             bias_mu = c(bias_mu_MM, bias_mu_PM, bias_mu_ML),
                             bias_sigma = c(bias_sigma_MM, bias_sigma_PM,
                                            bias_sigma_ML))

### 2.3.4 MSE calculation ----------------------------------------------------

MSE_mu_MM <- var(MM[,1])+bias_mu_MM^2
MSE_mu_PM <- var(PM[,1])+bias_mu_PM^2
MSE_mu_ML <- var(ML[,1])+bias_mu_ML^2

MSE_sigma_MM <- var(MM[,2])+bias_sigma_MM^2
MSE_sigma_PM <- var(PM[,2])+bias_sigma_PM^2
MSE_sigma_ML <- var(ML[,2])+bias_sigma_ML^2

results_normal <- cbind(results_normal,
                        MSE_mu = c(MSE_mu_MM, MSE_mu_PM, MSE_mu_ML),
                        MSE_sigma = c(MSE_sigma_MM, MSE_sigma_PM,
                                      MSE_sigma_ML))
results_normal

### 2.3.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 5000
step <- 3
size <- seq(size_min, size_max, step)

MM <- matrix(0,length(size),2,dimnames = list(c(),c("mu", "sigma")))
PM <- matrix(0,length(size),2,dimnames = list(c(),c("mu", "sigma")))
ML <- matrix(0,length(size),2,dimnames = list(c(),c("mu", "sgima")))
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rnorm(size[i], param[1], param[2])
  MM[i,] <- optim(c(4000,4000), MM_normal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  PM[i,] <- optim(c(4000,4000), PM_normal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  ML[i,] <- optim(c(4000,4000), ML_normal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  setTxtProgressBar(bar, i)
}

#### 2.3.5.1 mu parameter ----------------------------------------------------

MM[,1] <- MM[,1]-param[1]
PM[,1] <- PM[,1]-param[1]
ML[,1] <- ML[,1]-param[1]

plot(size, PM[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,1], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

#### 2.3.5.2 sigma parameter -------------------------------------------------

MM[,2] <- MM[,2]-param[2]
PM[,2] <- PM[,2]-param[2]
ML[,2] <- ML[,2]-param[2]

plot(size, PM[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,2], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

## 2.4 Log Normal ------------------------------------------------------------

### 2.4.1 Error functions ----------------------------------------------------

MM_lognormal<- function(param, sample){
  r1 <- (mean(log(sample)) - param[1])^2
  r2 <- (var(log(sample)) - param[2])^2
  r1 + r2
}

PM_lognormal <- function(param, sample){
  r1 <- (quantile(sample, 0.8)-qlnorm(0.8, param[1], param[2]))^2
  r2 <- (quantile(sample, 0.2)-qlnorm(0.2, param[1], param[2]))^2
  r1 + r2
}

ML_lognormal <- function(param, sample){
  -sum(dlnorm(sample, param[1], param[2], log = TRUE))
}

### 2.4.2 Parameters calculation ---------------------------------------------

n <- length(sev)
param <- c(8,1)
sample_DGP <- rlnorm(n, param[1], param[2])
mean(sample_DGP)
sd(sample_DGP)

optim(c(5,0.5), MM_lognormal, method = "L-BFGS-B", sample = sample_DGP)
optim(c(5,0.5), PM_lognormal, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))
optim(c(5,0.5), ML_lognormal, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))

### 2.4.3 Bias calculation ---------------------------------------------------

sim_size <- 500
MM <- matrix(0,sim_size,2,dimnames = list(c(),c("meanlog", "sdlog")))
PM <- matrix(0,sim_size,2,dimnames = list(c(),c("meanlog", "sdlog")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("meanlog", "sdlog")))

n <- length(sev)
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rlnorm(n, param[1], param[2])
  MM[i,] <- optim(c(5,0.5), MM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  PM[i,] <- optim(c(5,0.5), PM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(5,0.5), ML_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}
bias_meanlog_MM <- mean(MM[,1])-param[1]
bias_meanlog_PM <- mean(PM[,1])-param[1]
bias_meanlog_ML <- mean(ML[,1])-param[1]

bias_sdlog_MM <- mean(MM[,2])-param[2]
bias_sdlog_PM <- mean(PM[,2])-param[2]
bias_sdlog_ML <- mean(ML[,2])-param[2]

results_lognormal <- data.frame(method = c("MM", "PM", "ML"),
                             bias_meanlog = c(bias_meanlog_MM, bias_meanlog_PM,
                                              bias_meanlog_ML),
                             bias_sdlog = c(bias_sdlog_MM, bias_sdlog_PM,
                                            bias_sdlog_ML))

### 2.4.4 MSE calculation ----------------------------------------------------

MSE_meanlog_MM <- var(MM[,1])+bias_meanlog_MM^2
MSE_meanlog_PM <- var(PM[,1])+bias_meanlog_PM^2
MSE_meanlog_ML <- var(ML[,1])+bias_meanlog_ML^2

MSE_sdlog_MM <- var(MM[,2])+bias_sdlog_MM^2
MSE_sdlog_PM <- var(PM[,2])+bias_sdlog_PM^2
MSE_sdlog_ML <- var(ML[,2])+bias_sdlog_ML^2

results_lognormal <- cbind(results_lognormal,
                        MSE_meanlog = c(MSE_meanlog_MM, MSE_meanlog_PM,
                                        MSE_meanlog_ML),
                        MSE_sdlog = c(MSE_sdlog_MM, MSE_sdlog_PM,
                                      MSE_sdlog_ML))
results_lognormal

### 2.4.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 5000
step <- 3
size <- seq(size_min, size_max, step)

MM <- matrix(0,length(size),2,dimnames = list(c(),c("meanlog", "sdlog")))
PM <- matrix(0,length(size),2,dimnames = list(c(),c("meanlog", "sdlog")))
ML <- matrix(0,length(size),2,dimnames = list(c(),c("meanlog", "sdlog")))
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rlnorm(size[i], param[1], param[2])
  MM[i,] <- optim(c(5,0.5), MM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  PM[i,] <- optim(c(5,0.5), PM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, upper = c(20,2),
                  lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(5,0.5), ML_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}

#### 2.4.5.1 meanlog parameter -----------------------------------------------

MM[,1] <- MM[,1]-param[1]
PM[,1] <- PM[,1]-param[1]
ML[,1] <- ML[,1]-param[1]

plot(size, PM[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,1], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

#### 2.4.5.2 sdlog parameter -------------------------------------------------

MM[,2] <- MM[,2]-param[2]
PM[,2] <- PM[,2]-param[2]
ML[,2] <- ML[,2]-param[2]

plot(size, PM[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,2], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

## 2.5 Gamma -----------------------------------------------------------------

### 2.5.1 Error functions ----------------------------------------------------

# Due to convergence problems, the moments are directly calculated instead of
# using an error function as for the others. (Note the scale parameter theta
# is being used for the parameters estimation)

PM_gamma <- function(param, sample){
  r1 <- (quantile(sample, 0.8)-qgamma(0.8, param[1], scale = param[2]))^2
  r2 <- (quantile(sample, 0.2)-qgamma(0.2, param[1], scale = param[2]))^2
  r1 + r2
}

ML_gamma <- function(param, sample){
  -sum(dgamma(sample, param[1], scale = param[2], log = TRUE))
}

### 2.5.2 Parameters calculation ---------------------------------------------

n <- length(sev)
param <- c(5,0.5)
sample_DGP <- rgamma(n, param[1], scale = param[2])
mean(sample_DGP)
sd(sample_DGP)

c(mean(sample_DGP)^2/var(sample_DGP),var(sample_DGP)/mean(sample_DGP))
optim(c(4,1), PM_gamma, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))
optim(c(4,1), ML_gamma, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))

### 2.5.3 Bias calculation ---------------------------------------------------

sim_size <- 500
MM <- matrix(0,sim_size,2,dimnames = list(c(),c("alpha", "theta")))
PM <- matrix(0,sim_size,2,dimnames = list(c(),c("alpha", "theta")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("alpha", "theta")))

n <- length(sev)
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rgamma(n, param[1], scale = param[2])
  MM[i,] <- c(mean(sample_DGP)^2/var(sample_DGP),
              var(sample_DGP)/mean(sample_DGP))
  PM[i,] <- optim(c(4,1), PM_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(4,1), ML_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}
bias_alpha_MM <- mean(MM[,1])-param[1]
bias_alpha_PM <- mean(PM[,1])-param[1]
bias_alpha_ML <- mean(ML[,1])-param[1]

bias_theta_MM <- mean(MM[,2])-param[2]
bias_theta_PM <- mean(PM[,2])-param[2]
bias_theta_ML <- mean(ML[,2])-param[2]

results_gamma <- data.frame(method = c("MM", "PM", "ML"),
                                bias_alpha = c(bias_alpha_MM,
                                                 bias_alpha_PM,
                                                 bias_alpha_ML),
                                bias_sdlog = c(bias_theta_MM, bias_theta_PM,
                                               bias_theta_ML))

### 2.5.4 MSE calculation ----------------------------------------------------

MSE_alpha_MM <- var(MM[,1])+bias_alpha_MM^2
MSE_alpha_PM <- var(PM[,1])+bias_alpha_PM^2
MSE_alpha_ML <- var(ML[,1])+bias_alpha_ML^2

MSE_theta_MM <- var(MM[,2])+bias_theta_MM^2
MSE_theta_PM <- var(PM[,2])+bias_theta_PM^2
MSE_theta_ML <- var(ML[,2])+bias_theta_ML^2

results_gamma <- cbind(results_gamma,
                      MSE_alpha = c(MSE_alpha_MM, MSE_alpha_PM,
                                    MSE_alpha_ML),
                      MSE_theta = c(MSE_theta_MM, MSE_theta_PM,
                                    MSE_theta_ML))
results_gamma

### 2.5.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 5000
step <- 3
size <- seq(size_min, size_max, step)

MM <- matrix(0,length(size),2,dimnames = list(c(),c("alpha", "theta")))
PM <- matrix(0,length(size),2,dimnames = list(c(),c("alpha", "theta")))
ML <- matrix(0,length(size),2,dimnames = list(c(),c("alpha", "theta")))
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rgamma(size[i], param[1], scale = param[2])
  MM[i,] <- c(mean(sample_DGP)^2/var(sample_DGP),
              var(sample_DGP)/mean(sample_DGP))
  PM[i,] <- optim(c(4,1), PM_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(4,1), ML_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}

#### 2.5.5.1 alpha parameter -------------------------------------------------

MM[,1] <- MM[,1]-param[1]
PM[,1] <- PM[,1]-param[1]
ML[,1] <- ML[,1]-param[1]

plot(size, PM[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,1], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,1], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

#### 2.5.5.2 theta parameter -------------------------------------------------

MM[,2] <- MM[,2]-param[2]
PM[,2] <- PM[,2]-param[2]
ML[,2] <- ML[,2]-param[2]

plot(size, PM[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(2,1))
plot(size, MM[,2], type = "l")
abline(h = 0, col = "red")
plot(size, ML[,2], type = "l")
abline(h = 0, col = "red")
par(mfrow=c(1,1))

# 3. Distribution estimation -------------------------------------------------

## 3.1 Frecuency -------------------------------------------------------------

hist(frec)

### 3.1.1 As poisson ---------------------------------------------------------

param_pois <- optim(c(1), MM_pois, method = "L-BFGS-B", sample = frec)$par

#### 3.1.1.1 Graphical comparison --------------------------------------------

F_empiric_raw <- function(x, sample=sample_DGP){
  out <- c()
  for (i in 1:length(x)){
    out[i] <- mean(sample<=x[i])
  }
  out
}

dx <- 1E-5
range <- c(-1, -dx, dx)
for (i in 1:10) {
  range[2*i+2]=c(i-dx)
  range[2*i+3]=c(i+dx)
}
empiric <- F_empiric_raw(range,frec)
pois <- ppois(range, param_pois)

plot(range, empiric, type = "l", col = "blue")
lines(range, pois, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric", "Poisson Estimated"),
       col=c("blue", "red"), lty=1)

#### 3.1.1.2 Tests -----------------------------------------------------------

#install.packages("dgof")
library("dgof")

#install.packages("goftest")
library("goftest")

ks.test(frec,rpois(1000000,param_pois))

CVM_test_discrete <- function(sample, pFUN, param, range = c(0,10)){
  if (length(param)==1){
    cvm <- c(0)
    range <- seq(range[1], range[2])
    for (i in 1:length(range)){
      empiric <- mean(sample<=range[i])
      estimated <- pFUN(range[i],param)
      cvm <- cvm + (empiric-estimated)^2
    }
  }
  if (length(param)==2){
    cvm <- c(0)
    range <- seq(range[1], range[2])
    for (i in 1:length(range)){
      empiric <- mean(sample<=range[i])
      estimated <- pFUN(range[i],param[1], param[2])
      cvm <- cvm + (empiric-estimated)^2
    }
  }
  cvm*length(sample)
}

AD_test_discrete <- function(sample, pFUN, param, range = c(0,12)){
  if (length(param)==1){
    ad <- c(0)
    range <- seq(range[1], range[2])
    for (i in 1:length(range)){
      empiric <- mean(sample<=range[i])
      estimated <- pFUN(range[i],param)
      ad <- ad +
        (empiric-estimated)^2/(pFUN(range[1],param)*(1-pFUN(range[1],param)))
    }
  }
  if (length(param)==2){
    ad <- c(0)
    range <- seq(range[1], range[2])
    for (i in 1:length(range)){
      empiric <- mean(sample<=range[i])
      estimated <- pFUN(range[i],param[1], param[2])
      ad <- ad +
        (empiric-estimated)^2/(pFUN(range[1],param[1],
                                    param[2])*(1-pFUN(range[1],param[1],
                                                      param[2])))
    }
  }
  ad*length(sample)
}

cvm_pois <- c()
for (i in 1:5000){
  sample <- rpois(1116,param_pois)
  cvm_pois[i] <- CVM_test_discrete(sample, ppois, param_pois)
}
hist(cvm_pois)
quantile(cvm_pois, 0.95)
cvm_frec <- CVM_test_discrete(frec,ppois, param_pois)
mean(cvm_pois>cvm_frec)

ad_pois <- c()
for (i in 1:5000){
  sample <- rpois(1116,param_pois)
  ad_pois[i] <- AD_test_discrete(sample, ppois, param_pois)
}
hist(ad_pois)
quantile(ad_pois, 0.95)
ad_frec <- AD_test_discrete(frec,ppois, param_pois)
mean(ad_pois>ad_frec)

range <- seq(0,8,1)
expected <- dpois(range, param_pois)*length(frec)
obtained <- c()
for (i in 1:length(range)) obtained[i] = sum(frec==range[i])
chi_data_pois <- matrix(0,2,length(range))
chi_data_pois[1,] <- obtained
chi_data_pois[2,] <- expected
chisq.test(chi_data_pois)

### 3.1.2 As negative binomial -----------------------------------------------

param_nbinom <- optim(c(9,0.4), ML_nbinom, method = "L-BFGS-B", sample = frec,
                      upper = c(NA,0.9999), lower = c(0.00001,0.0001))$par

#### 3.1.2.1 Graphical comparison --------------------------------------------

dx <- 1E-5
range <- c(-1, -dx, dx)
for (i in 1:10) {
  range[2*i+2]=c(i-dx)
  range[2*i+3]=c(i+dx)
}
empiric <- F_empiric_raw(range,frec)
nbinom <- pnbinom(range, param_nbinom[1], param_nbinom[2])

plot(range, empiric, type = "l", col = "blue")
lines(range, nbinom, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric", "Negative Binomial Estimated"),
       col=c("blue", "red"), lty=1)

#### 3.1.2.2 Tests -----------------------------------------------------------

ks.test(frec,rnbinom(10000000, param_nbinom[1], param_nbinom[2]))

cvm_nbinom <- c()
for (i in 1:5000){
  sample <- rnbinom(1116,param_nbinom[1], param_nbinom[2])
  cvm_nbinom[i] <- CVM_test_discrete(sample, pnbinom, param_nbinom)
}
hist(cvm_nbinom)
quantile(cvm_nbinom, 0.95)
cvm_nbinom_frec <- CVM_test_discrete(frec,pnbinom, param_nbinom)
mean(cvm_nbinom>cvm_nbinom_frec)

ad_nbinom <- c()
for (i in 1:5000){
  sample <- rpois(1116,param_pois)
  ad_nbinom[i] <- AD_test_discrete(sample, pnbinom, param_nbinom)
}
hist(ad_nbinom)
quantile(ad_nbinom, 0.95)
ad_frec <- AD_test_discrete(frec,pnbinom, param_nbinom)
mean(ad_nbinom>ad_frec)

range <- seq(0,8,1)
expected <- dnbinom(range, param_nbinom[1], param_nbinom[2])*length(frec)
obtained <- c()
for (i in 1:length(range)) obtained[i] = sum(frec==range[i])
chi_data_nbinom <- matrix(0,2,length(range))
chi_data_nbinom[1,] <- obtained
chi_data_nbinom[2,] <- expected
chisq.test(chi_data_nbinom)

## 3.2 Severity --------------------------------------------------------------

hist(sev, breaks = 20)

### 3.2.1 As Normal ----------------------------------------------------------

param_normal <- optim(c(4000,4000), ML_normal,
                      method = "L-BFGS-B", sample = sev)$par

#### 3.2.1.1 Graphical comparison --------------------------------------------

range <- seq(-10000, 40000)

empiric <- F_empiric_raw(range,sev)
norm <- pnorm(range, param_normal[1],param_normal[2])

plot(range, empiric, type = "l", col = "blue")
lines(range, norm, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric", "Normal"),
       col=c("blue", "red"), lty=1)

# it can not be normal since negative severity values are not possible

### 3.2.2 As Log Normal ------------------------------------------------------

param_lnorm <- optim(c(8,1.3), ML_lognormal, method = "L-BFGS-B", sample = sev,
               lower = c(0.001,0.001))$par

#### 3.2.2.1 Graphical comparison --------------------------------------------

range <- seq(-5000, 30000)

empiric <- F_empiric_raw(range,sev)
lnorm <- plnorm(range, param_lnorm[1],param_lnorm[2])

plot(range, empiric, type = "l", col = "blue")
lines(range, lnorm, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric", "Normal"),
       col=c("blue", "red"), lty=1)

#### 3.2.2.2 Tests -----------------------------------------------------------

ks.test(sev,rlnorm(10000000, param_lnorm[1], param_lnorm[2]))

sim <- 1000
cvm_results <- matrix(0,sim,2)
ad_results <- matrix(0,sim,2)
for (i in 1:sim){
  out <- cvm.test(sev, plnorm, meanlog = param_lnorm[1],
                  sdlog = param_lnorm[2], estimated = TRUE)
  cvm_results[i,] <- c(out$statistic, out$p.value)
  out <- ad.test(sev, plnorm, meanlog = param_lnorm[1],
                 sdlog = param_lnorm[2], estimated = TRUE)
  ad_results[i,] <- c(out$statistic, out$p.value)
  
}
mean(cvm_results[,1])
mean(cvm_results[,2])
mean(ad_results[,1])
mean(ad_results[,2])

### 3.2.3 As Gamma -----------------------------------------------------------

param_gamma <- optim(c(0.8,100), ML_gamma, method = "L-BFGS-B", sample = sev,
                    lower = c(0.001,0.001))$par

#### 3.2.3.1 Graphical comparison --------------------------------------------

range <- seq(-5000, 30000)

empiric <- F_empiric_raw(range,sev)
gamma <- pgamma(range, param_gamma[1], scale = param_gamma[2])

plot(range, empiric, type = "l", col = "blue")
lines(range, gamma, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric", "Gamma"),
       col=c("blue", "red"), lty=1)

#### 3.2.3.2 Tests -----------------------------------------------------------

ks.test(sev,rgamma(10000000, param_gamma[1],scale = param_gamma[2]))

sim <- 1000
cvm_results <- matrix(0,sim,2)
ad_results <- matrix(0,sim,2)
for (i in 1:sim){
  out <- cvm.test(sev, pgamma, shape = param_gamma[1],
                             scale = param_gamma[2], estimated = TRUE)
  cvm_results[i,] <- c(out$statistic, out$p.value)
  out <- ad.test(sev, pgamma, shape = param_gamma[1],
                           scale = param_gamma[2], estimated = TRUE)
  ad_results[i,] <- c(out$statistic, out$p.value)

}
mean(cvm_results[,1])
mean(cvm_results[,2])
mean(ad_results[,1])
mean(ad_results[,2])

# 4. Aggregated model for non-life -------------------------------------------

policies <- 25234

sim <- 10000
cost_nl <- c()
bar <- txtProgressBar(0,sim,style=3)
for (i in 1:sim){
  n <- sum(rpois(policies, param_pois))
  cost_nl[i] <- sum(rgamma(n,param_gamma[1], shape = param_gamma[2]))
  setTxtProgressBar(bar, i)
}
hist(cost_nl, breaks = 25)
est_cost_nl <- mean(cost_nl)
VaR995_nl <- quantile(cost_nl,0.995)
TVaR995_nl <- mean(cost_nl[cost_nl>=VaR995_nl])

# 5. Life model --------------------------------------------------------------

mort$qx <- as.double(mort$qx)

sample_size <- 10000
deceased_sample <- matrix(0, length(mort$Edad),sample_size,
                      dimnames = list(mort$Edad))

for (i in 1:length(mort[,1])){
  deceased_sample[i,] <- rbinom(sample_size,size = mort$ni[i],prob = mort$qx[i])
}

# we check if the dead rate is as expected
mean(deceased_sample[2,])/mort$ni[2]
mort$qx[2]

sim <- 10000
deceased <- c()
bar <- txtProgressBar(0,sim,style=3)
for (i in 1:sim){
  total <- c()
  for (j in 1:length(mort[,1])){
    total[j] <- sample(deceased_sample[j,],1)
  }
  deceased[i] <- sum(total)
  setTxtProgressBar(bar, i)
}
amount <- 1E6
cost_l <- deceased*amount
hist(cost_l, breaks = 25)
VaR995_deceased <- quantile(deceased, 0.995)
TVaR995_deceased <- mean(deceased[deceased>VaR995_deceased])
est_cost_l <- mean(cost_l)
VaR995_l <- VaR995_deceased*amount
TVaR995_l <- TVaR995_deceased*amount

# 6. Life and Non-life aggregated model --------------------------------------

sim <- 10000

# assuming independence between variables
cost_total <- c()
for (i in 1:sim){
  cost_total[i] <- sample(cost_l,1)+sample(cost_nl,1)
}

hist(cost_total, breaks = 25)
est_cost_total <- mean(cost_total)
VaR995_total <- quantile(cost_total,0.995)
TVaR995_total <- mean(cost_total[cost_total>=VaR995_total])