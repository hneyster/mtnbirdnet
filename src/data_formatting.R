library(here)
library(magrittr)

dat <- read.csv(here("raw/seymour.csv"))

# turn this csv file into an array with dimensions sites X reps x species

alt <- dat$elevation 

sites <- unique(dat$location) %>% sort()
nsite <- length(sites) 
spec <- unique(dat$Scientific.Name) %>% as.character()
nspec <- length(spec)
rep<- c("a","b")
nrep <- length(rep)
yc<-array(0,dim = c(nsite, nrep, nspec),dimnames = list("site"=sites, "rep"=rep, "spec"=spec))

for (k in 1:nspec){
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (l in 1:nrow(cvld)) {
        yc[i,j,k] <- ifelse(cvld$pc.number[l]==sites[i] && cvld$visit.x[l]==rep[j] && 
                            cvld$species[l]==spec[k], 
                          cvld$number[l], 
                          yc[i,j,k])
      }
    }
  }
}
ycc<-ifelse(is.na(yc), 0,yc) ## taking out the NAs to zeros 

nsite <- dat$Location %>% unique %>% length()
nsp <- dat$Scientific.Name %>% unique %>% length()
nrep <- 2

y <-array(0,dim = c(nsite, nrep, nspec),dimnames = list("site"=sites, "rep"=rep, "spec"=spec))
