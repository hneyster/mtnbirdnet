# R code for post processing of stan model of elevation responses 
# Written by Harold Eyster github @hneyster 

library(shinystan)
library(magrittr)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(bayesplot)
library(rstan)
library(bayestestR)
library(cowplot)
library(here)
ci_fxn <- function(x) ci(x,method="HDI", ci=0.50) #gives 89% CrI (using HDI)
iter = 4000
load(here('fmt/pcdat.Rdata'))
spnames<- dimnames(pcdat$y)[3]$spec

#load("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211115/")
csvs <- c("/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211115/output_1.csv",
          "/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211115/output_2.csv",
          "/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211115/output_3.csv",
          "/home/harold/Dropbox/gitfiles/seymour/fit_seymour_20211115/output_4.csv")
#fit <- as_cmdstan_fit(csvs, format = "draws_matrix")
fit <- as_cmdstan_fit(csvs)

#fitcsv <-read_cmdstan_csv(csvs)

draws <- fit$draws(c("M","a","b","c","d"), format = "draws_matrix")
draws_array <- array(dim=c(iter, pcdat$S,5 ), data = draws)

mdraw <- fit$draws("M",format ='draws_matrix') %>% as.data.frame()
adraw <- fit$draws("a",format ='draws_matrix') %>% as.data.frame()
bdraw <- fit$draws("b",format ='draws_matrix') %>% as.data.frame()
cdraw <- fit$draws("c",format ='draws_matrix') %>% as.data.frame()
ddraw <- fit$draws("d",format ='draws_matrix') %>% as.data.frame()

msum <- fit$summary("M", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
asum <- fit$summary("a", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
bsum <- fit$summary("b", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
csum <- fit$summary("c", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
dsum <- fit$summary("d", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()

pc <- fit$summary("N", 'mean', q25=~quantile2(., probs = c(0.25, 0.75))) %>% as.data.frame()
#dpc <- fit$draws("N", format = 'draws_matrix') %>% as.data.frame()
#apc <- array(dim = c(iter, pcdat$S, pcdat$R), data = dpc)
matpc <- matrix(ncol = pcdat$S, nrow= pcdat$R, data = pc$mean)
lpc <- matrix(ncol = pcdat$S, nrow= pcdat$R, data = pc$q25)
hpc <- array(dim = c(pcdat$R,pcdat$S, 3), data = pc$q75)



pc_ci <- NA
for (i in 1:5){
  tmp<- do.call("rbind", apply(draws_array[,,i], FUN=ci_fxn, MARGIN = 2)) %>% data.frame() %>% 
    cbind('median' = apply(draws_array[,,i], FUN =median, MARGIN = 2))  %>% cbind( "sp"=spnames) %>% cbind('var' = paste(i))
  draws_ci <- rbind(draws_ci, tmp)
} 


X11()
pm <- mcmc_intervals(draws,pars = vars(starts_with("M["))) + scale_y_discrete(labels = (dimnames(pcdat$y)[3]$spec)) + xlab("M")
pa <- mcmc_intervals(draws,pars = vars(starts_with("a["))) + scale_y_discrete(labels = (dimnames(pcdat$y)[3]$spec)) + xlab("a")+ yaxis_text(on=FALSE)
pb <- mcmc_intervals(draws,pars = vars(starts_with("b["))) + scale_y_discrete(labels = (dimnames(pcdat$y)[3]$spec)) + xlab("b") + yaxis_text(on=FALSE)
pc <- mcmc_intervals(draws,pars = vars(starts_with("c["))) + scale_y_discrete(labels = (dimnames(pcdat$y)[3]$spec)) + xlab("c") + yaxis_text(on=FALSE)
pd <- mcmc_intervals(draws,pars = vars(starts_with("d["))) + scale_y_discrete(labels = (dimnames(pcdat$y)[3]$spec)) + xlab("d") + yaxis_text(on=FALSE)
paramplot <- plot_grid(pm,pa,pb,pc,pd, ncol = 5, rel_widths = c(2,1,1,1,1))
pdf(file = here("out/paramplot_1121.pdf"),width = 12, height = 10)
paramplot
dev.off()
poisson_log_rng(M[s]*(1/(1+exp(a[s]+b[s]*elevation_std[i])))*(1/(1+exp(c[s]+d[s]*elevation_std[i]))))
x <- pcdat$elevation_std

resp ={ function(x, sp, metric) exp(msum[sp,metric]*
                             (1/(1+exp(asum[sp,metric]+bsum[sp,metric]*x)))*
                             (1/(1+exp(csum[sp,metric]+dsum[sp,metric]*x))))}

resp_pred <- data.frame("elevation"=pcdat$elevation, "elevation_std"=pcdat$elevation_std)

resp_plot <- vector('list', length(spnames))
for (sp in 1:length(spnames)){
  resp_plot[[sp]]  <- local({
    sp <- sp
    resp_pred$y <- lapply(resp_pred$elevation_std, FUN = function(x) resp(x,sp, metric=2))
    resp_pred$l <- lapply(resp_pred$elevation_std, FUN = function(x) resp(x,sp, metric=3))
    resp_pred$h <- lapply(resp_pred$elevation_std, FUN = function(x) resp(x,sp, metric=4))
    
   plotsp<-  ggplot(resp_pred) +
     geom_ribbon(aes(x=as.double(elevation), y=as.double(y), ymin =as.double(l), 
   ymax = as.double(h)), color = 'gray',fill='gray')+ geom_line(aes(x=as.double(elevation), y=as.double(y)), color = 'black')+
     ggtitle(paste(spnames[sp])) + ylab('pred. abund.') + xlab('elevation (m)') + theme_bw() + theme(plot.title = element_text(size = 10))
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