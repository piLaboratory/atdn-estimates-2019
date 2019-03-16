data {
    real<lower=0> a; // fraction of community that has been sampled
    real<lower=0> lN; // log of the total number of individuals in the whole system
    int<lower=1> s; // number of species in the sample
    int<lower=1> y[s]; // abundances of species in the sample
    }

    parameters {
      real mu;
      real tau;
      //real<lower=0> lam[s];
      real llam[s]; //log of abundances
      real ephi;
    }

    transformed parameters {
        real<lower=0> sig;
        real<lower=0> phi;
        sig = exp(tau);
        phi = exp(ephi);
    }

    model {
      mu ~ normal(15, 30);
      tau ~ normal(1, 2);
      ephi ~ normal(0, 5);
      llam ~ normal(mu, sig);
      for(i in 1:s){
	//llam[i] ~ normal(mu, sig);
	y[i] ~ neg_binomial_2_log(llam[i]+log(a), phi);//first argument is log(mean)
	target += -log1m((phi/(exp(llam[i]+log(a))+phi))^phi);// Adjusts log-likelihodd for zero-truncation in  negative binomial (http://discourse.mc-stan.org/t/error-with-truncation-on-negative-binomial-distribution-log-alternative-parameterization/1051/2)
        }
    }

   generated quantities {
     real S1; // derived quantity (Estimated number of species from the lognormal that has been sampled by negative binomial)
     S1 = exp(lN - (mu + sig^2/2) );
}

