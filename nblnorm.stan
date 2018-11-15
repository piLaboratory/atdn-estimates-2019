data {
    real<lower=0> a; // fraction of community that has been sampled
    real<lower=0> lN; // log of the total number of individuals in the whole system
    int<lower=1> s; // number of species in the sample
    int<lower=0> y[s]; // abundances of species in the sample
    }

    parameters {
        real mu;
        real tau;
        real<lower=0> lam[s];
        real ephi;
    }

    transformed parameters {
        real<lower=0> sig;
        real<lower=0> phi;
        sig = exp(tau);
        phi = exp(ephi);
    }

    model {
        mu ~ normal(10, 30);
        tau ~ normal(1, 3);
        ephi ~ normal(0, 5);
        for(i in 1:s){
          lam[i] ~ lognormal(mu, sig);
          y[i] ~ neg_binomial_2(lam[i]*a, phi);
        }
    }

   generated quantities {
        real S1; // derived quantity (Estimated number of species from the lognormal that has been sampled by negative binomial)
        S1 = exp(lN - (mu + sig^2/2) );
}

