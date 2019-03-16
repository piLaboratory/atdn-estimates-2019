/* Zero-truncated Beta-binomial distribution as a explicit hierarchical model
a.k.a a compound beta-binomial distribution truncated at zero.
 */
data {
  int<lower=1> Nsites; // total number of sites
  int<lower=1> Sobs; // observed number of species
  int<lower=0> y[Sobs]; // number sites each species has been recorded
}
parameters {
  real<lower=0,upper=1> p[Sobs]; // probability of occurrence per plot for each species recorded
  real<lower=0,upper=1> phi; // transformed parameters of beta, see Stan manual p 286
  real<lower=0.1> lambda;	
}

transformed parameters{
  real alpha;
  real beta;
  alpha = lambda * phi;
  beta = lambda*(1-phi);
}

model {
  phi ~ beta(1,1); //prior recommended in Stan manual p.286
  lambda ~ pareto(0.1, 1.5); //prior recommended in Stan manual p.286
  for(i in 1:Sobs){
    p[i] ~ beta(alpha, beta);
    y[i] ~ binomial(Nsites, p[i])T[1,]; // zero-truncated binomial
  }
}

generated quantities {
  // total number of species from an untruncated beta-binomial
  real Sest;
  Sest = Sobs / (1 - beta_binomial_cdf(1 , Nsites, alpha, beta));
}

