library(rstan)
library(here)
library(magrittr)
stdz <- function(x){
  y <- (x-mean(x))/(2*sd(x))
} # for stdzing by 2 std deviations. See Gelman 2008 


dat <- read.csv(here("clean/seymour.csv")) # the cleaned data
# turn this csv file into an array with dimensions sites X reps x species

# removing point counts only visited once:  
locations_by_checklist <- dat[,c("Submission.ID","Location.ID")]
locations_by_checklist_un <- dat[!duplicated(locations_by_checklist),] #data frame with only unique combinations of location and checklist 
dat <- dat[dat$Location.ID %in% names(which(table(locations_by_checklist_un$Location.ID) > 1)), ] #removing the locations that were only visited once 

# selecting the first point count: 
locations_by_checklist_g1 <- dat[,c("Submission.ID","Location.ID")] # locations by checklists of ones visited more than once 
locations_by_checklist_g1_un <- locations_by_checklist_g1[!duplicated(locations_by_checklist_g1),] #data frame with only unique combinations of location and checklist 

dat$rep<-0
rep1 <- !duplicated(locations_by_checklist_g1_un$Location.ID) #first visit
dat[dat$Submission.ID %in% locations_by_checklist_g1_un$Submission.ID[rep1],"rep"] <- 1
rep2 <- !duplicated(locations_by_checklist_g1_un$Location.ID,fromLast = TRUE) #last visit
dat[dat$Submission.ID %in% locations_by_checklist_g1_un$Submission.ID[rep2],"rep"] <-2

dat <- dat[!dat$rep==0,] #removing additional visits that occurred after the first two visits. 


sub <- unique(dat$Submission.ID)
nsub <- length(sub)
site <- unique(dat$Location.ID) 
nsite <- length(site) 
spec <- unique(dat$Common.Name) %>% as.character()
nspec <- length(spec)
rep<- unique(dat$rep)
nrep <- length(rep)
y<-array(0,dim = c(nsite, nrep, nspec),dimnames = list("site"=site, "rep"=rep, "spec"=spec))

for (k in 1:nspec){
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (l in 1:nrow(dat)) {
        y[i,j,k] <- ifelse(dat$Location.ID[l]==site[i] && dat$rep[l]==rep[j] && 
                            dat$Common.Name[l]==spec[k], 
                          dat$Count[l], 
                          y[i,j,k])
      }
    }
  }
}

y<-ifelse(is.na(y), 0,y) ## NA = 0

elevation <- dat[!duplicated(dat$Location.ID), "elevation"] #getting elevation for each site 
elevation_std <- stdz(elevation)

pcdat<-list(y=y, elevation = elevation, elevation_std = elevation_std, R=nsite, J=nrep, S=nspec)

save(pcdat, file = here("fmt/pcdat.Rdata"))
HETH <-  list(y=y[,,"Hermit Thrush"], elevation = elevation, elevation_std = elevation_std, R=nsite, J=nrep)
SWTH <- list(y=y[,,"Swainson's Thrush"], elevation = elevation, elevation_std = elevation_std, R=nsite, J=nrep)
save(HETH, file = here("fmt/HETH.Rdata"))
save(SWTH, file = here("fmt/SWTH.Rdata"))
obs <- y
y<- obs[,,"Hermit Thrush"]
K <- 10
R = nsite
J = nrep
stan_rdump(list =  c('y','elevation', 'elevation_std','R','J','K'), file = here("fmt/HETH.Rdump"))

y <- obs[,,"Swainson's Thrush"]
stan_rdump(list =  c('y','elevation', 'elevation_std','R','J','K'), file = here("fmt/SWTH.Rdump"))
