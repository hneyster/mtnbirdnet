# R code for post processing of stan model of elevation responses 
# Written by Harold Eyster github @hneyster 

library(magrittr)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(bayesplot)
library(rstan)
library(bayestestR)
library(cowplot)
library(here)
ci_fxn <- function(x) ci(x,method="HDI", ci=0.89) #gives 89% CrI (using HDI)
iter = 4000
spfile <- "HETH"
rundate <- '20220322'
load(here('fmt/SWTH.Rdata'))
load(here('fmt/HETH.Rdata'))
R <- SWTH$R
x <- SWTH$elevation_std
xreal <- SWTH$elevation
S <- 1


singleplot <- function(spfile,data,rundate='20220322') {
max <- apply(data$y, FUN =max, MARGIN =1)
csvs <- c(
    paste("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_",rundate,"/",spfile,"_1.csv", sep = ""),
    paste("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_",rundate,"/",spfile,"_2.csv", sep = ""),
    paste("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_",rundate,"/",spfile,"_3.csv", sep = ""),
    paste("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_",rundate,"/",spfile,"_4.csv", sep = ""))
         
#fit <- as_cmdstan_fit(csvs, format = "draws_matrix")
fit <- as_cmdstan_fit(csvs)

#fitcsv <-read_cmdstan_csv(csvs)

draws <- fit$draws(c("M","a","b","c","d","p"), format = "draws_matrix")
draws_array <- array(dim=c(iter, 1,5 ), data = draws)

m <- fit$draws("M",format ='draws_matrix') %>% unlist() %>% as.double()
a <- fit$draws("a",format ='draws_matrix')  %>% unlist() %>% as.double()
b <- fit$draws("b",format ='draws_matrix')  %>% unlist() %>% as.double()
c <- fit$draws("c",format ='draws_matrix')  %>% unlist() %>% as.double()
d <- fit$draws("d",format ='draws_matrix')  %>% unlist() %>% as.double()

pred <- matrix(nrow=iter, ncol = R)
abund <- function(x, m.vec=m, a.vec=a,b.vec=b,c.vec=c,d.vec=d)  (m.vec*(1/(1+exp(a.vec+b.vec*x)))* (1/(1+exp(c.vec-d.vec*x))))
for (i in 1:R){
  pred[,i] <- abund(x[i])
}
pred <- as.data.frame(pred)
names(pred) <-xreal
sums <- ci_fxn(pred) %>% as.data.frame()
sums$mean <- colMeans(pred)
predplot<- ggplot(sums, aes(x= xreal, y = mean, ymin = CI_low, ymax = CI_high)) +geom_ribbon(fill = 'gray') + geom_line( color = 'black')+ geom_point(aes(y=max))+
    ylab(paste(spfile,' predicted abundance (89% CrI) \ndots = data',sep="")) + xlab('elevation (m)') + theme_bw() + theme(plot.title = element_text(size = 10))

return(predplot)
}
param <- mcmc_intervals(draws,pars = c('M','a','b','c','d','p'))

HETHplot <- singleplot(spfile = "HETH",data = HETH)
SWTHplot <- singleplot("SWTH",SWTH)
pdf(file = here('out/HETH.SWTH.pdf'))
plot_grid(SWTHplot, HETHplot, ncol=1)
dev.off()
#############################################
############################################
########### RXV  ############################
############################################
############################################



# mcmc_intervals(pred) + scale_y_discrete(labels = HETH$elevation) + coord_flip()
msum <- fit$summary("M", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
asum <- fit$summary("a", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
bsum <- fit$summary("b", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
csum <- fit$summary("c", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
dsum <- fit$summary("d", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()

# pc <- fit$summary("N", 'mean', q25=~quantile2(., probs = c(.055, 0.945))) %>% as.data.frame()
#dpc <- fit$draws("N", format = 'draws_matrix') %>% as.data.frame()
#apc <- array(dim = c(iter, pcdat$S, pcdat$R), data = dpc)
#matpc <- matrix(ncol = S, nrow= R, data = pc$mean)
#lpc <- matrix(ncol = S, nrow= R, data = pc$q25)
#hpc <- array(dim = c(R,S, 3), data = pc$q75)


ggplot(pc, aes(x=SWTH$elevation, y=mean, ymin=q5.5, ymax=q94.5)) + geom_ribbon(fill='gray70') + geom_line() + geom_point(y=max)
plot(SWTH$elevation,pc$mean)
plot(SWTH$elevation, max)



# pdf(file = here("out/paramplot_1121.pdf"),width = 12, height = 10)
#paramplot
#dev.off()
#poisson_log_rng(M[s]*(1/(1+exp(a[s]+b[s]*elevation_std[i])))*(1/(1+exp(c[s]+d[s]*elevation_std[i]))))

resp ={ function(x, metric) exp(msum[metric]*
                             (1/(1+exp(asum[metric]+bsum[metric]*x)))*
                             (1/(1+exp(csum[metric]-dsum[metric]*x))))}

resp_pred <- data.frame("elevation"=SWTH$elevation, "elevation_std"=SWTH$elevation_std)

#resp_plot <- vector('list', length(spnames))
#for (sp in 1:length(spnames)){
#  resp_plot[[sp]]  <- local({
#    sp <- sp
    resp_pred$y <- lapply(resp_pred$elevation_std, FUN = function(x) resp(x, metric=2))
    resp_pred$l <- lapply(resp_pred$elevation_std, FUN = function(x) resp(x, metric=3))
resp_pred$h <- lapply(resp_pred$elevation_std, FUN = function(x) resp(x, metric=4))

resp_pred <- apply(resp_pred, FUN = function(x) (as.double(unlist(x))), MARGIN = 2) %>% as.data.frame()
   plotsp<-  ggplot(resp_pred) + geom_ribbon(aes(x=elevation, y =y, ymin = l, ymax = h), color = 'gray',fill='gray')+ geom_line(aes(x=elevation, y=(y)), color = 'black')+
      ylab('pred. abund.') + xlab('elevation (m)') + theme_bw() + theme(plot.title = element_text(size = 10))
    print(plotsp)
  })
}
p1<- plot_grid(plotlist = resp_plot[c(1:47)], ncol = 5)#, labels = c(spnames[1:47]), label_y = 0.47, label_x = .3, label_fontface = 'plain')
pdf(file = here("out/elevation_response_1121.pdf"),width = 9, height = 15)
p1
dev.off()


for (sp in 1:pcdat$S){
  # abunp[,sp] <- exp(mean(mdraw[,s]*(1/(1+exp(adraw[,s]+bdraw[,s]*pcdat$elevation_std)))*(1/(1+exp(cdraw[,s]+ddraw[,s]*elevation_std)))))
  
  plot(pcdat$elevation, matpc[,sp],
       main = paste(dimnames(pcdat$y)[3]$spec[sp]),
       xlab = "elevation",
       xlim = c(200, 1400),
       ylim = c(0,5))
}
mcmc_dens(as.data.frame(apc[,1,]))

#############################################
############################################
########### RXV  ############################
############################################
############################################

#load("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211022.Rdata")

#launch_shinystan(fit_seymour_20211022)
draws <- rstan::extract(fit_seymour_20211022)
launch_shinystan(fit)

pdf(file = here("out/elevation_response_1121.pdf"),width = 9, height = 15)
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
