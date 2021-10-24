# R code for post processing of stan model of elevation responses 
# Written by Harold Eyster github @hneyster 

library(shinystan)
library(magrittr)
library(here)

load("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211022.Rdata")

#launch_shinystan(fit_seymour_20211022)
draws <- rstan::extract(fit_seymour_20211022)

pdf(file = here("out/elevation_response.pdf"),width = 9, height = 15)
par(mar=c(2,2,2,2))
par(mfrow = c(10,5))
abunp <- matrix(ncol = pcdat$S, nrow=pcdat$R, NA) 
for (sp in 1:pcdat$S){
  abunp[,sp] <- exp(mean(draws$alpha[,sp]) + mean(draws$beta1[,sp])*pcdat$elevation_std + mean(draws$beta2[,sp])*pcdat$elevation_std^2)
  
  plot(pcdat$elevation, abunp[,sp],
       main = paste(dimnames(pcdat$y)[3]$spec[sp]),
       xlab = "elevation",
       xlim = c(200, 1400),
       ylim = c(0,5))
}
dev.off()

par(mfrow = c(1,1))
plot(pcdat$elevation, rowSums(abunp),
     main = "Number of individuals",
     xlab = "elevation",
     xlim = c(200, 1400),
     ylim = c(0,40))


abunp_bin <- (abunp > .9)*1L
pdf(file = here("out/species_richness.pdf"))
plot(pcdat$elevation, rowSums(abunp_bin),
     ylab = "Species richness",
     xlab = "elevation",
     xlim = c(200, 1400),
     ylim = c(0,40))
dev.off()
# 
# pdf(file = here("out/elevation_response_abs.pdf"),width = 9, height = 15)
# par(mar=c(2,2,2,2))
# par(mfrow = c(10,5))
# abun <- matrix(ncol = pcdat$S, nrow=pcdat$R, NA) 
# for (sp in 1:pcdat$S){
#   abun[,sp] <- (mean(draws$alpha[,sp]) + mean(draws$beta1[,sp])*pcdat$elevation_std + mean(draws$beta2[,sp])*pcdat$elevation_std^2)
#   
#   plot(pcdat$elevation, abun[,sp],
#        main = paste(dimnames(pcdat$y)[3]$spec[sp]),
#        xlab = "elevation",
#        xlim = c(200, 1400),
#        ylim = c(-4,2))
# }
# dev.off()
#sp1 <- plot(pcdat$elevation,(mean(draws$alpha[,sp]) + mean(draws$beta1[,sp])*pcdat$elevation_std + mean(draws$beta2[,sp])*pcdat$elevation_std^2),ylab = paste(dimnames(pcdat$y)[3]$spec[sp]), xlab="elevation")
#plot(pcdat$elevation, sp1)