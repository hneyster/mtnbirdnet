## Function to simulate bird communities, for testing in Bayesian multispecies
## abundance models. 

## Written by Harold Eyster haroldeyster@gmail.com 
## Github @hneyster
## March 2022
library(rstan)
library(here)
## First singlespecies:
stdz <- function(x){
  y <- (x-mean(x))/(2*sd(x))
}
R = 400
x_meters <- sample(1:1400, R)
x <- stdz(x_meters)
J=2
m = 4
a = -1
b = 6
c = -2
d= 2
p <- .7
abund <- function(x)  (m*(1/(1+exp(a+b*x)))* (1/(1+exp(c-d*x))))
y <- abund(x)
#plot(x_meters,y)
z <- NA
for (i in 1:length(x)){z[i] <- rpois(n=1, lambda=y[i])}
b <- matrix(nrow=R, ncol = J, NA)
for (i in 1:R) (b[i,] <- rbinom(n=2, size = z[i], prob = p))

y <- b
K<-15
elevation_std <- x
elevation <- x_meters
stan_rdump(list =  c('y','elevation', 'elevation_std','R','J','K'), file = here("fmt/sim_1sp.Rdump"))

################# R X V #########################

simbird <- function (type = c("det/nondet", "counts"), nfarms = 10, nsite = 100, nrep = 2, 
          nspec = 60, mean.lambda = 2, sig.loglam = 1, mu.beta.loglam = 1, 
          sig.beta.loglam = 1, sig.farm = 1, mean.p = 0.25, sig.lp = 1, mu.beta.lp = 0, 
          sig.beta.lp = 0, show.plot = TRUE, seed= NULL) 
{
  if (FALSE) 
    x <- NULL
  type <- match.arg(type)
  if (show.plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    oldAsk <- devAskNewPage(ask = TRUE)
    on.exit(devAskNewPage(oldAsk), add = TRUE)
  }
  set.seed(seed)
    if (type == "counts") {
      nsites.per.farm <- nsite/nfarms
      y.all <- y.obs <- p <- array(NA, c(nsite, nrep, nspec))
      dimnames(y.all) <- dimnames(y.obs) <- dimnames(p) <- list(paste("site", 
             1:nsite, sep = ""), paste("rep", 1:nrep, sep = ""), 
             paste("sp", 1:nspec, sep = ""))
      N <- lambda <- matrix(NA, nsite, nspec)
      dimnames(N) <- dimnames(lambda) <- list(paste("site", 
          1:nsite, sep = ""), paste("sp", 1:nspec, sep = ""))
      detected.at.all <- rep(NA, nspec)
      habitat <- rep(0:1,nsites.per.farm/2) #making either zero or one. For three habitats: 
                                          #sample(1:3,10, replace=TRUE)
      wind <- matrix(rnorm(nsite * nrep), ncol = nrep)
      mu.loglam <- log(mean.lambda)
      mu.lp <- ifelse(mean.p == "1", 500, qlogis(mean.p))
      farm_log_lam<-rnorm(nfarms, mu.loglam,sig.farm)
      beta0<-array(NA, dim = c(nspec, nfarms))
        for (f in 1:nfarms) {
          beta0[,f] <-rnorm(nspec, farm_log_lam[f], sig.loglam)
        }
      # beta0<-array(NA, dim = nsites, 1)
      # counter<-0
      #  for (f in 1:nfarms) {
      #     beta0[counter+1:nsites.per.farm] <-rnorm(nsites.per.farm, farm_log_lam[f], sig.loglam)
      #      counter<-nsites.per.farm*f
      #   }
      # 
      # beta0 <- rnorm(nspec, mu.loglam, sig.loglam)
      beta1 <- rnorm(nspec, mu.beta.loglam, sig.beta.loglam)
      alpha0 <- rnorm(nspec, mu.lp, sig.lp)
      alpha1 <- rnorm(nspec, mu.beta.lp, sig.beta.lp)
      for (k in 1:nspec) {
        counter<-0
        for (f in 1:nfarms) {
        lambda[(counter+1):(nsites.per.farm*f), k] <- exp(beta0[k,f] + beta1[k] * habitat)
        counter<-nsites.per.farm*f
          for (j in 1:nrep) {
            p[, j, k] <- plogis(alpha0[k] + alpha1[k] * wind[, j])
          }
        }
      }
      for (k in 1:nspec) {
        N[, k] <- rpois(nsite, lambda[, k])
      }
      tmp <- apply(N, 2, sum)
      occurring.in.sample <- as.numeric(tmp > 0)
      for (k in 1:nspec) {
        for (i in 1:nsite) {
          for (j in 1:nrep) {
            y.all[i, j, k] <- rbinom(1, N[i, k], p[i, j, k])
          }
        }
        detected.at.all[k] <- if (any(y.all[, , k] > 0)) 
          TRUE
        else FALSE
      }
      y.obs <- y.all[, , detected.at.all]
      detected.at.site <- apply(y.obs > 0, c(1, 3), any)
      ymax.obs <- apply(y.all, c(1, 3), max)
      Ntotal.fs <- sum(occurring.in.sample)
      Ntotal.obs <- sum(detected.at.all)
      
      if (show.plot) {
        par(mfrow = c(1, 2), mar = c(5, 5, 5, 3), cex.axis = 1.3, 
            cex.lab = 1.3)
        curve(exp(beta0[1] + beta1[1] * x), -2, 2, main = "Species-specific (black) and community (red) \n response of lambda to habitat", 
              xlab = "Habitat", ylab = "Expected abundance (lambda)")
        for (k in 1:nspec) {
          curve(exp(beta0[k] + beta1[k] * x), -2, 2, add = TRUE)
        }
        curve(exp(mu.loglam + mu.beta.loglam * x), -2, 2, 
              col = "red", lwd = 3, add = TRUE)
        curve(plogis(alpha0[1] + alpha1[1] * x), -2, 2, main = "Species-specific (black) and community (red) \n response of detection to wind", 
              xlab = "Wind", ylab = "Detection probability (p)", 
              ylim = c(0, 1))
        for (k in 2:nspec) {
          curve(plogis(alpha0[k] + alpha1[k] * x), -2, 
                2, add = TRUE)
        }
        curve(plogis(mu.lp + mu.beta.lp * x), -2, 2, col = "red", 
              lwd = 3, add = TRUE)
        par(mfrow = c(2, 2), mar = c(5, 6, 4, 2), cex.axis = 1.3, 
            cex.lab = 1.3)
        mapPalette <- colorRampPalette(c("yellow", "orange", 
                                         "red"))
        image(x = 1:nspec, y = 1:nsite, z = log10(t(N) + 
                                                    1), col = mapPalette(100), main = paste("True log(abundance) (log10(N)) matrix\n (finite-sample N species =", 
                                                                                            sum(occurring.in.sample), ")"), frame = TRUE, 
              xlim = c(0, nspec + 1), zlim = c(0, log10(max(N))), 
              xlab = "Species", ylab = "Sites")
        image(x = 1:nspec, y = 1:nsite, z = log10(t(ymax.obs) + 
                                                    1), col = mapPalette(100), main = paste("Observed maximum counts (log10 + 1)"), 
              xlim = c(0, nspec + 1), frame = TRUE, xlab = "Species", 
              ylab = "Sites", zlim = c(0, log10(max(N))))
        ratio <- ymax.obs/N
        ratio[ratio == "NaN"] <- 1
        image(x = 1:nspec, y = 1:nsite, z = t(ratio), col = mapPalette(100), 
              main = paste("Ratio of max count to true abundance (N)"), 
              xlim = c(0, nspec + 1), frame = TRUE, xlab = "Species", 
              ylab = "Sites", zlim = c(0, 1))
        lims <- c(0, log10(max(N + 1)))
        plot(log(N), log(ymax.obs), xlab = "True abundance (log10(N+1))", 
             ylab = "Observed max count \n(log10(max+1))", 
             xlim = lims, ylim = lims, main = "Observed vs. true N (log10 scale)")
        abline(0, 1)
      }
  return(list(type = type, nsite = nsite, nrep = nrep, 
              nspec = nspec, nfarms=nfarms, mean.lambda = mean.lambda, mu.loglam = mu.loglam, 
              sig.loglam = sig.loglam, mu.beta.loglam = mu.beta.loglam, 
              sig.beta.loglam = sig.beta.loglam, farm_log_lam = farm_log_lam, sig.farm = sig.farm, mean.p = mean.p, 
              mu.lp = mu.lp, sig.lp = sig.lp, mu.beta.lp = mu.beta.lp, 
              sig.beta.lp = sig.beta.lp, lambda = lambda, p = p, beta0=beta0, 
              N = N, y.all = y.all, y.obs = y.obs, ymax.obs = ymax.obs, 
              Ntotal.fs = Ntotal.fs, Ntotal.obs = Ntotal.obs))
}
}
