
// simulating the prior-predictive distribution for the "...spfarm.stan" program
// super helpful post: https://discourse.mc-stan.org/t/choosing-weakly-informative-priors-for-population-level-effects-in-a-poisson-glmm/18008/11
data {
  int<lower=0> R;       // Number of sites
  int<lower=0> J;       // Number of temporal replicates
  int<lower=0> S;       // Number of species 
  int<lower=0> y[R, J, S]; // Counts
  real elevation_std[R];          // elevation of each point count
}
 
generated quantities{
 // Params
vector[S] alpha;
vector[S] beta1;
vector[S] beta2;

real mu_alpha;
real mu_beta1;
real mu_beta2;
real <lower=0> sig_alpha;
real <lower=0> sig_beta1;
real <lower=0> sig_beta2;

real  mu_logit_p;
real <lower=0> sig_logit_p;
vector[S] logit_p; // detection probability by species 

real logN[R,S];
  int N[R,S]; 

// priors 
  mu_alpha = normal_rng(0,1);
  sig_alpha = fabs(normal_rng(0,2));
  mu_beta1 = normal_rng(0,1);
  sig_beta1 = fabs(normal_rng(0,2));
  mu_beta2 = normal_rng(-1,1);
  sig_beta2 = fabs(normal_rng(0,2));
 
mu_logit_p  =  normal_rng(0, 2); 
sig_logit_p =  fabs(normal_rng(1, 3));


  for (s in 1:S){
    alpha[s] =  normal_rng(mu_alpha, sig_alpha); 
    beta1[s] =  normal_rng(mu_beta1, sig_beta1); 
    beta2[s] =  normal_rng(mu_beta1, sig_beta2); 
    logit_p[s] =   logistic_rng(mu_logit_p, sig_logit_p);
  }

  for (i in 1:R){
    for (s in 1:S){
        logN[i, s] = alpha[s] + beta1[s]*elevation_std[i] + beta2[s]*(elevation_std[i])^2; 
    N[i,s] = poisson_log_rng(logN[i,s]);
    }
  }
}
