library(rstan)
library(mvtnorm)
library(xtable)

rstan_options(auto_write = TRUE) # save compiled model

load("array_allcont")

set.seed(123)
initdays <- 14
N <- 9
CONT <- 4
N1 <- 2
N2 <- 4
N3 <- 2
N4 <- 1
Tl <- 290 - initdays
C <- 5


p1dirich <- rep(0.5, 5) + apply(array_allcont[1:initdays,,],3,sum)
y <- array_allcont[-(1:initdays),,]

par(mfrow=c(3,3))
for (i in 1:9) {
  matplot(y[,i,], main=i)
}

data = list(N = N, CONT = CONT, N1 = N1, N2 = N2, N3 = N3, N4 = N4, T = Tl, C = C, y = y, p1dirich = p1dirich)
stan_fit = stan("model_sensitivity_1_p1.stan",
                data = data,
                iter = 5000,
                warmup = 2500,
                refresh = 5,
                chains = 4, cores =4, control = list(max_treedepth = 17))

print(stan_fit, pars=c('alpha', 's0', 'p1', 'Sigma_NA', 'Omega_NA','Sigma_EU', 'Omega_EU','Sigma_AS', 'Omega_AS'),
      digits = 5)
