//  binomial mixture model
// written by Harold Eyster github.com/hneyster
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
   // following Huisman et al. (1993) and Oksanen and Minchen (2002), the HOF model is: y = M11 + exp(a + bx)11 + exp(c-dx)
      
//    vector <lower=0> [S]  sig_farm;     // sd of pop sizes across farms 
    
//    vector[H] mu_habitat;
//    vector<lower=0>[H] sig_habitat;
//    matrix [H,S] habitat_tilde; // centered habitat weights 
    
//    matrix [F,S] farm_tilde;
    vector [S] logit_p_tilde;  // centered detection variation by species 

    real mu_alpha;
    real <lower=0> sig_alpha; 
    vector[S] alpha_tilde; //centered abundance by species
    
    real mu_beta1;
    real <lower=0> sig_beta1; 
    vector[S] beta1_tilde; //centered abundance by species
    
    real mu_beta2;
    real <lower=0> sig_beta2; 
    vector[S] beta2_tilde; //centered abundance by species
}
transformed parameters{
//  matrix[H,S] habitat; // habitat sensitivity by species 
  vector[S] logit_p; // detection probability by species 
//  matrix[F,S] farm; // farm sensitivity by species
//  for (i in 1:H){
//    for (j in 1:S){
 //   habitat[i,j] = mu_habitat[i]+sig_habitat[i]*habitat_tilde[i,j];
  //implies habitat ~ normal(mu_habitat, sig_habitat
//    }
//  }
  logit_p = mu_logit_p+sig_logit_p*logit_p_tilde;
  // implies logit_p ~ logistic(mu_logit_p, sig_logit_p)
//    for (i in 1:F){
//     for (j in 1:S){
//      farm[i,j] = sig_farm[j]*farm_tilde[i,j];
//      //implies farm ~ normal(0, sig_farm) 
//       }
//     }
 vector[S] alpha;
 alpha = mu_alpha+sig_alpha*alpha_tilde;
 
  vector[S] beta1;
 beta1 = mu_beta1+sig_beta1*beta1_tilde;

   vector[S] beta2;
 beta2 = mu_alpha+sig_beta2*beta2_tilde;
}
model {
  // Priors
  mu_alpha ~ normal(0,5);
  sig_alpha ~ normal(0,2);
  mu_beta1 ~ normal(0,5);
  sig_beta1 ~ normal(0,2);
  mu_beta2 ~ normal(0,5);
  sig_beta2 ~ normal(0,2);
  // habitat is: 1= corn, 2=soy, 3=hay, 4=prarie, 5=woods, 6=young polyculture, 7= old polyculture 
//  mu_habitat[1] ~normal(-2, 1); //corn
//  mu_habitat[2] ~normal(-2, 1); //soy
//  mu_habitat[3] ~normal(-2, 1); //hay

//  mu_habitat[4] ~ normal(-1,1); // prarie
//  mu_habitat[5] ~ normal(-1,1); // wood 
  
//  mu_habitat[6] ~ normal(-1,1); // young poly
//  mu_habitat[7] ~ normal(-1,1); // mature poly 
  
//  sig_habitat ~ normal(0,2);

//  sig_farm ~ normal(0, 1); // we expect small differences between farms 
  mu_logit_p ~ normal(0,5); // prior for the mean detection probability, in logistic space (=0.5 in real scale). 
  sig_logit_p ~ normal(1,3); // scale parameter for logistic dist.  
  
  logit_p_tilde ~ logistic(0,1); // implies logit_p ~ logistic (mu_logit_p, sig_logit_p) because  logistic(mu, sig) = logistic(0,1)*sig + mu. see https://mc-stan.org/users/documentation/case-studies/mle-params.html

for (i in 1:S){
  alpha[i] ~ std_normal();
  beta1[i] ~ std_normal();
  beta2[i] ~ std_normal();
}

//  for (i in 1:F){
//   for (j in 1:S){
//    farm_tilde[i,j] ~ std_normal(); 
//      //implies farm ~ normal(0, sig_farm) 
//    }
//  }
//  for (i in 1:H){
//    for (j in 1:S){
//  habitat_tilde[i,j] ~ std_normal(); //implies habitat ~ normal(mu_habitat, sig_habitat)
//    }
// }
  
  // Likelihood
  for (i in 1:R) {
    for (s in 1:S){
      vector[K - max_y[i,s] + 1] lp; //it's the product, lambda*p that's calculated
      
      for (j in 1:(K - max_y[i,s] + 1))
        lp[j] = poisson_log_lpmf(max_y[i,s] + j - 1 | alpha[s] + beta1[s]*elevation_std[i] + beta2[s]*(elevation_std[i])^2) //lpmf bc poisson is discrete. 
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
       N[i,s] = poisson_log_rng(alpha[s] + beta1[s]*elevation_std[i] + beta2[s]*(elevation_std[i])^2);
//       //PPC[i,s] = binomial_rng (N[i,s], logit_p[s] ) ;
     }
   }
 }
