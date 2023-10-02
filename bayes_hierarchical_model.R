## bayesian hierarchical model

library(rstan)
library(mvtnorm)
library(ggplot2)
library(xtable)
rstan_options(auto_write = TRUE) # save compiled model

#### sampling for the model through the stan file
load("array_allcont.rda")
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
p1dirich <- rep(1,5) + apply(array_allcont[1:initdays,,],3,sum)
y <- array_allcont[-(1:initdays),,]

data = list(N = N, CONT = CONT, N1 = N1, N2 = N2, N3 = N3, N4 = N4, T = Tl, C = C, y = y, p1dirich = p1dirich)
stan_fit = stan("bayes_hierarchical_model.stan",
                data = data,
                iter = 5000,
                warmup = 2500,
                refresh = 5,
                chains = 4, cores =4, control = list(max_treedepth = 17))


#### print the posterior parameters 
print(stan_fit, pars=c('alpha', 's0', 'p1', 'Sigma_NA', 'Omega_NA',
                       'Sigma_EU', 'Omega_EU','Sigma_AS', 'Omega_AS'))
a <- summary(stan_fit, pars=c('alpha', 's0', 'p1'))$summary
b <- summary(stan_fit, pars=c('Sigma_NA','Sigma_EU', 'Sigma_AS'))$summary
c <- summary(stan_fit, pars=c('Omega_NA','Omega_EU', 'Omega_AS'))$summary

xtable(a[,c(1,4,8)] ,digits=c(4,4,4,4))
xtable(b[,c(1,4,8)] ,digits=c(4,4,4,4))
xtable(c[,c(1,4,8)] ,digits=c(4,4,4,4))

#### use the samples to plot the posterior means and 95% credible intervals for each country
param <- extract(stan_fit)
p1 <- extract(stan_fit, pars=c('p1'))[["p1"]]
alpha <- extract(stan_fit, pars=c('alpha'))[["alpha"]]
eps_NA <- extract(stan_fit, pars=c('eps_NA'))[["eps_NA"]]
eps_EU <- extract(stan_fit, pars=c('eps_EU'))[["eps_EU"]]
eps_AS <- extract(stan_fit, pars=c('eps_AS'))[["eps_AS"]]
eps_AU <- extract(stan_fit, pars=c('eps_AU'))[["eps_AU"]]

samp <- 10000
p <- array(dim=c(samp,276,9,5))
p_avg <- array(dim=c(samp,276,9,5))
y_sim <- array(dim=c(samp,276,9,5))

for (i in 1:N){
  sapply(1:samp, function(k) p[k,1,i,1:5] <<- p1[k, 1:5])
}


for (k in 1:samp){
  alpha_k <- alpha[k,1:4]
  lr_j <- numeric(4)
  
  for (t in 2:Tl){
    for (i in 1:N){
      if (i<=N1){
        
        for (j in 2:C) {
          lr_j[j-1] = log(p[k,t-1,i,j]/p[k,t-1,i,1]) + alpha_k[j-1] + eps_NA[k,t-1,j-1,i]
        }
        
        #print(t)
        
        pexp = exp(lr_j)
        p[k,t,i,1] = 1/(sum(pexp)+1)
        

        for (j in 2:C) {
          p[k,t,i,j] = p[k,t,i,1] * pexp[j-1]
        }
      }
      
      if (i>N1 && i<=N1+N2){
        
        for (j in 2:C) {
          lr_j[j-1] = log(p[k,t-1,i,j]/p[k,t-1,i,1]) + alpha_k[j-1] + eps_EU[k,t-1,j-1,i-N1]
        }
        pexp = exp(lr_j)
        p[k,t,i,1] = 1/(sum(pexp)+1)
        
        for (j in 2:C) {
          p[k,t,i,j] = p[k,t,i,1] * pexp[j-1]
        }
      }
      
      if (i>N1+N2 && i<=N1+N2+N3){
        
        for (j in 2:C) {
          lr_j[j-1] = log(p[k,t-1,i,j]/p[k,t-1,i,1]) + alpha_k[j-1] + eps_AS[k,t-1,j-1,i-N1-N2]
        }
        pexp = exp(lr_j)
        p[k,t,i,1] = 1/(sum(pexp)+1)
        
        for (j in 2:C) {
          p[k,t,i,j] = p[k,t,i,1] * pexp[j-1]
        }
      }
      
      if (i>N1+N2+N3){
        
        for (j in 2:C) {
          lr_j[j-1] = log(p[k,t-1,i,j]/p[k,t-1,i,1]) + alpha_k[j-1] + eps_AU[k,t-1,j-1]
        }
        pexp = exp(lr_j)
        p[k,t,i,1] = 1/(sum(pexp)+1)
        
        for (j in 2:C) {
          p[k,t,i,j] = p[k,t,i,1] * pexp[j-1]
        }
      }
    }
  }
}

set.seed(100)

for (k in 1:samp){
  for (i in 1:N){
    for (t in 2:Tl){
      size <- sum(y[t,i,1:C])
      
      if (t < 15){
        for (j in 1:C){
          p_avg[k,t,i,j] <- mean(p[k, 1:t-1,i,j])
        }
      } else{
        for (j in 1:C){
          p_avg[k,t,i,j] <- mean(p[k, (t-14):(t-1),i,j])
        }
      }
      prob <- p_avg[k,t,i,1:C]
      y_sim[k,t,i,1:C] <- rmultinom(n=1,size=size, prob = prob)
      
    }
  }
}

praw <- array(dim=c(276,9,5))

for (i in 1:N){
  for (t in 1:Tl){
    praw[t,i,1:C] <- y[t,i,1:C]/sum(y[t,i,1:C])
  }
}

#### posterior means and 95% credible intervals of p
pmean <- array(dim=c(276,9,5))
plow <- array(dim=c(276,9,5))
phigh <- array(dim=c(276,9,5))

for (i in 1:N){
  for (t in 2:Tl){
    for (j in 1:C){
      pmean[t,i,j] <- mean(p_avg[1:samp,t,i,j])
      plow[t,i,j] <- quantile(p_avg[1:samp,t,i,j], 0.025)
      phigh[t,i,j] <- quantile(p_avg[1:samp,t,i,j], 0.975)
    }
  }
}

#### compute posterior means and 95% credible intervals of y
ymean <- array(dim=c(276,9,5))
ylow <- array(dim=c(276,9,5))
yhigh <- array(dim=c(276,9,5))

for (i in 1:N){
  for (t in 2:Tl){
    for (j in 1:C){
      ymean[t,i,j] <- mean(y_sim[1:samp,t,i,j])
      ylow[t,i,j] <- quantile(y_sim[1:samp,t,i,j], 0.025)
      yhigh[t,i,j] <- quantile(y_sim[1:samp,t,i,j], 0.975)
    }
  }
}


#### plot Figure 4
clist <- c("United States (US)", "Canada (CA)", "United Kingdom (UK)", "Netherlands (NL)", 
           "France (FR)", "Spain (SP)", "China (CN)", "India (IN)", "Australia (AU)")
clist <- factor(clist, levels=c("United States (US)", "Canada (CA)", "United Kingdom (UK)", "Netherlands (NL)", 
                                "France (FR)", "Spain (SP)", "China (CN)", "India (IN)", "Australia (AU)"))
raw <- numeric()
est <- numeric()
low <- numeric()
high <- numeric()
for (i in 1:9){
  raw <- c(raw, c(praw[,i,1], praw[,i,2], praw[,i,3], praw[,i,4], praw[,i,5]))
  est <- c(est, c(pmean[,i,1], pmean[,i,2], pmean[,i,3], pmean[,i,4], pmean[,i,5]))
  low <- c(low, c(plow[,i,1], plow[,i,2], plow[,i,3], plow[,i,4], plow[,i,5]))
  high <- c(high, c(phigh[,i,1], phigh[,i,2], phigh[,i,3], phigh[,i,4], phigh[,i,5]))
}

df <- data.frame(Date = rep(rep(seq(as.Date("2020-01-07"), as.Date("2020-10-08"), by="days"), 5),9),
                 raw = raw,
                 est = est,
                 low = low,
                 high = high,
                 category = rep(c("Cluster I","Cluster II", "Cluster III", "Cluster IV", "Cluster V"), 
                                each = 276),
                 country = rep(clist, each = 276*5))

p <- ggplot(df, aes(x=Date, y = raw, color = category)) + geom_point() + 
  geom_line(aes(y=est))+ 
  geom_ribbon(aes(ymin = low, ymax = high, fill = category), alpha = 0.3) + 
  xlab("Dates") + ylab("proportions") + ylim(c(0,1)) + facet_wrap(~country) +
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())



pdf(file="plot_posterior_cluster.pdf", width = 9, height = 8)  
p
dev.off()


