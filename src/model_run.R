library(rstan)
library(here)
options(mc.cores = parallel::detectCores())
load(here("fmt/pcdat.Rdata"))

fit_seymour_20211115<-stan(here("src/binmix_elevation.stan"), 
              data = c(pcdat, K=40), #max count in the data is 10 
             chains = 1, iter = 500, save_warmup=FALSE,
             seed = 1, control=list(max_treedepth=15)) #init_r=1, control=list(adapt_delta=0.9))


 save(fit_seymour_20211115, file="/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211115.Rdata")

# fit_seymour_20211022<-stan(here("src/binmix_elevation.stan"), 
#                             data = c(pcdat, K=40), #max count in the data is 10 
#                             chains = 1, iter = 500, save_warmup=FALSE,
#                             seed = 1) #init_r=1, control=list(adapt_delta=0.9))
#  
#  
#save(fit_seymour_20211022, file="/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211022.Rdata")
 
 