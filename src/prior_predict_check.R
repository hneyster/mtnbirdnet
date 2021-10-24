library(magrittr)
library(rstan)
library(here)
load(here("fmt/pcdat.Rdata"))
prior_pred <- stan("prior_predic.stan", data = c(pcdat), algorithm="Fixed_param", seed=1, iter = 100, chains = 1)
list_of_draws <- rstan::extract(prior_pred)
ab <- list_of_draws$N
meanab <- apply(MARGIN = c(2:3), FUN=mean, ab)
summary(prior_pred, pars = "mu_logit_p")$summary
summary(prior_pred, pars = "mu_logit_p")
stan_hist(prior_pred, pars="mu_logit_p")
stan_hist(prior_pred, pars = "N[1,1]")
library(bayesplot)
sigp<- list_of_draws$sig_logit_p 
logit_p <- list_of_draws$logit_p[1]
mu_p <- list_of_draws$mu_logit_p
mu_p_trans <- plogis(mu_p)
p <- rnorm(2000, mu_p, sigp) # this has way too much weight around 0 and 1
sim_logit_p<- rlogis(2000, 0,.5)
rlogis(10000, 2,.5) %>% plogis() %>% hist()
 # if the shape paramter is > 1, then you get bimodality
 # sign of the shape parameter doesn't matter 
runif(1000, -.5,.5) %>% hist() # to select the mu_logit_p
rnorm(1000, 0,.5) %>% hist() # or if we want to avoid uniform 
runif(1000, .05, 1) # for the scale parameter, sig_logit_p 
rnorm(1000, .5, 0.25) %>% hist()

# testing
plogis(list_of_draws$logit_p[6,]) %>% hist()
list_of_draws$logit_p %>% dim()
hist(p)
hist(qlogis(
hist(plogis(p))
mcmc_hist(prior_pred, pars = "mu_logit_p")
, 
          transformations = list("mu_logit_p" = "exp"))
prior_pred <- stan("prior_predic.stan", data = c(pcdat_all), algorithm="Fixed_param", seed=1, iter = 1000, chains = 4)
summary(prior_pred, pars = "mu_logit_p")$summary
summary(prior_pred, pars = "mu_logit_p")
stan_hist(prior_pred, pars="mu_logit_p")
library(bayesplot)
list_of_draws <- rstan::extract(prior_pred)
sigp<- list_of_draws$sig_logit_p 
logit_p <- list_of_draws$logit_p[1]
mu_p <- list_of_draws$mu_logit_p
mu_p_trans <- plogis(mu_p)
p <- rnorm(2000, mu_p, sigp) # this has way too much weight around 0 and 1
sim_logit_p<- rlogis(2000, 0,.5)
rlogis(10000, 2,.5) %>% plogis() %>% hist()
 # if the shape paramter is > 1, then you get bimodality
 # sign of the shape parameter doesn't matter 
runif(1000, -.5,.5) %>% hist() # to select the mu_logit_p
rnorm(1000, 0,.5) %>% hist() # or if we want to avoid uniform 
runif(1000, .05, 1) # for the scale parameter, sig_logit_p 
rnorm(1000, .5, 0.25) %>% hist()

# testing
plogis(list_of_draws$logit_p[,78]) %>% hist()
list_of_draws$logit_p %>% dim()
hist(p)
hist(qlogis(
hist(plogis(p))


## testing hab + farm priors
muhab<- list_of_draws$mu_habitat
muhab %>% hist()
hab <- list_of_draws$habitat
farm <- list_of_draws$farm
hist(rpois(2000, (hab[,1,1])))
max (hab[,1,1])
max(hab) #20.7944

logN <- list_of_draws$logN # this is log(lambda) in the poisson model
exp(logN) %>% hist() # this is lambda
mcmc_hist(prior_pred, pars = "logN", transformations = list(log N = "exp"))
max (logN) # this is 29 -- way too big
max(pcdat_all$y) # 38
log(50) # we shouldn't expect things to get much bigger than this.
hist(pcdat_all$y, ylim=c(0,100))
hist(rpois(10000,exp(0)))
hist(logN)
hist(exp(logN), ylim= c(0,1000))
hist(rnorm(10000, 0, 2) + (rnorm(10000, 0,2) + rnorm(10000, 0,1))) # basically equivalent
rpois(10000,exp(rnorm(10000, 0, 2) + (rnorm(10000, 0,2) + rnorm(10000, 0,1)))) %>% hist

(exp(rnorm(10000, 0, 1) + (rnorm(10000, 0,1) + rnorm(10000, 0,1)))) %>% hist(ylim=c(0,100))
curve( dlnorm( x , -1 , 1+1 ) , from=0 , to=2 , n=200 )
curve( dlnorm( x , 1 , 2 ) , from=0 , to=100 , n=200 )
# mean = exp(μ + σ2 /2),
exp(0 + 2^2/2)
exp(0 + 1/2)
for (i in 1:7){
    mean(pcdat_all$y[pcdat_all$hab_index==i,,]) %>% print()
}
exp(1+4/2)
exp(1+1^2/2)

hist(list_of_draws$N, xlim=c(0,100)
list_of_draws$N %>% summary()
curve(dlogis(x, .5, .5), from = 0, to = 1, n =200)
rlogis(10000, .5, .5) %>% summary()
(((rlogis(10000, 0, 1))*.5 + .5)) %>% summary() ## proof that logistic(mu, sig) = logistic(0,1)*sig + mu 
