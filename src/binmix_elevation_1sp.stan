//  binomial mixture model
// written by Harold Eyster github.com/hneyster
// Now for just one species
// Mar 22, 2022
data {
  int<lower=0> R;       // Number of sites
  int<lower=0> J;       // Number of temporal replicates
  int<lower=0> y[R, J]; // Counts
  int<lower=0> K;       // Upper bound of population size
  real elevation_std[R];          // elevation of each point count
  
}

transformed data {
  int<lower=0> max_y[R]; // the the max number of each sp observed at each spatial replication.
  
  for (i in 1:R) {
      max_y[i] = max(y[i,]); // 
  }
}

parameters {
    real<lower=0, upper=1> p; 
    real<lower=0> M;
    real a;
    real b;
    real c;
    real d;
    

}

model {
  // Priors
  M ~ normal(0,5);
  a ~ normal(0,5);
  b ~ normal(0,5);
  c ~ normal(0,5);
  d ~ normal(0,5);

  p ~ normal(.5,.5); // prior for the detection probability
 

  // Likelihood
   // following Huisman et al. (1993) and Oksanen and Minchen (2002), the HOF model is: y = M*1(/(1 + exp(a + bx))*1/(1 + exp(c-dx)))

  for (i in 1:R) {
      vector[K - max_y[i] + 1] lp; //it's the product, lambda*p that's calculated
      
      for (j in 1:(K - max_y[i] + 1))
        lp[j] = poisson_lpmf(max_y[i] + j - 1 | (M*(1/(1+exp(a+b*elevation_std[i])))*(1/(1+exp(c-d*elevation_std[i]))))) //lpmf bc poisson is discrete. 
      // max_y[i]+j-1 is n_i
      //This reports the probability that n_i takes this value given different samples of lambda 
      + binomial_lpmf(y[i,] | max_y[i] + j - 1, p); // implicitly sums across T (y[i] is vectorized)
      target += log_sum_exp(lp);

  }
}

generated quantities{
  int N[R];
  // real pred[R];
 for (i in 1:R){
       N[i] = poisson_rng(M*(1/(1+exp(a+b*elevation_std[i])))*(1/(1+exp(c-d*elevation_std[i]))));
//       //PPC[i,s] = binomial_rng (N[i], logit_p ) ;

   
     }
 }
