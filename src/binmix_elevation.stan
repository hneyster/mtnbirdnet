//  binomial mixture model
// written by Harold Eyster github.com/hneyster
// adding non-gaussian responses
data {
  int<lower=0> R;       // Number of sites
  int<lower=0> J;       // Number of temporal replicates
  int<lower=0> S;       // Number of species 
  int<lower=0> y[R, J, S]; // Counts
  int<lower=0> K;       // Upper bound of population size
  real elevation_std[R];          // elevation of each point count
  
}

transformed data {
  int<lower=0> max_y[R,S]; // the the max number of each sp observed at each spatial replication.
  
  for (i in 1:R) {
    for (j in 1:S)
      max_y[i,j] = max(y[i,,j]); // 
  }
}

parameters {
    real mu_logit_p;       // 
    real <lower=0> sig_logit_p; // 
    vector [S] logit_p_tilde;  // centered detection variation by species 
    
    real mu_M;
    real <lower=0> sig_M; 
    vector[S] M_tilde; //centered abundance by species

    real mu_a;
    real <lower=0> sig_a; 
    vector[S] a_tilde; //centered abundance by species
    
    real mu_b;
    real <lower=0> sig_b; 
    vector[S] b_tilde; //centered abundance by species
    
    real mu_c;
    real <lower=0> sig_c; 
    vector[S] c_tilde; //centered abundance by species
    
    real mu_d;
    real <lower=0> sig_d; 
    vector[S] d_tilde; //centered abundance by species

}
transformed parameters{
  vector[S] logit_p; // detection probability by species 

  logit_p = mu_logit_p+sig_logit_p*logit_p_tilde;

 vector[S] M;
 M = mu_M+sig_M*M_tilde;
 
vector[S] a;
 a = mu_a+sig_a*a_tilde;

   vector[S] b;
 b = mu_b+sig_b*b_tilde;
 
 vector[S] c;
 c = mu_c+sig_c*c_tilde;
 
 vector[S] d;
 d = mu_d+sig_d*d_tilde;
}
model {
  // Priors
  mu_M ~ normal(0,5);
  sig_M ~ normal(0,2);
  mu_a ~ normal(0,5);
  sig_a ~ normal(0,2);
  mu_b ~ normal(0,5);
  sig_b ~ normal(0,2);
  mu_c ~ normal(0,5);
  sig_c ~ normal(0,2);
  mu_d ~ normal(0,5);
  sig_d ~ normal(0,2);

  mu_logit_p ~ normal(0,5); // prior for the mean detection probability, in logistic space (=0.5 in real scale). 
  sig_logit_p ~ normal(1,3); // scale parameter for logistic dist.  
  
  logit_p_tilde ~ logistic(0,1); // implies logit_p ~ logistic (mu_logit_p, sig_logit_p) because  logistic(mu, sig) = logistic(0,1)*sig + mu. see https://mc-stan.org/users/documentation/case-studies/mle-params.html

for (s in 1:S){
  M[s] ~ std_normal();
  a[s] ~ std_normal();
  b[s] ~ std_normal();
  c[s] ~ std_normal();
  d[s] ~ std_normal();
}

  // Likelihood
   // following Huisman et al. (1993) and Oksanen and Minchen (2002), the HOF model is: y = M*1(/(1 + exp(a + bx))*1/(1 + exp(c-dx))

  for (i in 1:R) {
    for (s in 1:S){
      vector[K - max_y[i,s] + 1] lp; //it's the product, lambda*p that's calculated
      
      for (j in 1:(K - max_y[i,s] + 1))
        lp[j] = poisson_log_lpmf(max_y[i,s] + j - 1 | (M[s]*(1/(1+exp(a[s]+b[s]*elevation_std[i])))*(1/(1+exp(c[s]-d[s]*elevation_std[i]))))) //lpmf bc poisson is discrete. 
      // max_y[i]+j-1 is n_i
      //This reports the probability that n_i takes this value given different samples of lambda 
      + binomial_logit_lpmf(y[i,,s] | max_y[i,s] + j - 1, logit_p[s]); // implicitly sums across T (y[i] is vectorized)
      target += log_sum_exp(lp);
    }
  }
}

generated quantities{
  int N[R,S];
//  // int PPC[R,S];
// 
 for (i in 1:R){
    for (s in 1:S){
       N[i,s] = poisson_log_rng(M[s]*(1/(1+exp(a[s]+b[s]*elevation_std[i])))*(1/(1+exp(c[s]-d[s]*elevation_std[i]))));
//       //PPC[i,s] = binomial_rng (N[i,s], logit_p[s] ) ;
     }
   }
 }
