// Not implemented yet. Expected number of species is a parameter to estimate but still integer
functions {

  /**
   * (from http://discourse.mc-stan.org/t/poisson-binomial-distribution-any-existing-stan-implementation/4220/7)
   * Return the Poisson-binomial log probability mass for the specified
   * count y and vector of probabilities theta.  The result is the log
   * probability of y successes in N = size(theta) independent
   * trials with probabilities of success theta[1], ..., theta[N].
   *
   * See:  https://en.wikipedia.org/wiki/Poisson_binomial_distribution
   *
   * @param y number of successes
   * @param theta vector of probabilities
   * @return Poisson-binomial log probability mass
   */
  real poisson_binomial_lpmf(int y, vector theta) {
    int N = rows(theta);
    matrix[N + 1, N + 1] alpha;
    vector[N] log_theta = log(theta);
    vector[N] log1m_theta = log1m(theta);

    if (y < 0 || y > N)
      reject("poisson binomial variate y out of range, y = ", y,
             " N = ", N);
    for (n in 1:N)
      if (theta[n] < 0 || theta[n] > 1)
        reject("poisson binomial parameter out of range,",
               " theta[", n, "] =", theta[n]);

    if (N == 0)
      return y == 0 ? 0 : negative_infinity();

    // dynamic programming with cells
    // alpha[n + 1, tot + 1] = log prob of tot successes in first n trials
    alpha[1, 1] = 0;
    for (n in 1:N) {
      // tot = 0
      alpha[n + 1, 1] = alpha[n, 1] + log1m_theta[n];

      // 0 < tot < n
      for (tot in 1:min(y, n - 1))
        alpha[n + 1, tot + 1]
            = log_sum_exp(alpha[n, tot] + log_theta[n],
                          alpha[n, tot  + 1] + log1m_theta[n]);

      // tot = n
      if (n > y) continue;
      alpha[n + 1, n + 1] = alpha[n, n] + log_theta[n];
    }
    return alpha[N + 1, y + 1];
  }

}
data {
  int<lower=1> Nsites; // total number of sites
  int<lower=1> Sobs; // observed number of species
  int<lower=0> y[Sobs]; // number sites each species has been recorded
  int<lower=Sobs> Smax; //Upper limit of number of species to test in brute force seek of true number of species (see derived quantities)
}
parameters {
  real<lower=0,upper=1> p[Sobs]; // probability of occurrence per plot for each species recorded
  real ls1; // Parameters of the beta distribution of probability of occurrences
  real ls2;	
}
transformed parameters {
  real<lower=0> s1;
  real<lower=0> s2;
  s1 = exp(ls1);
  s2 = exp(ls2);
}
model {
  ls1 ~ normal(-1, 100);
  ls2 ~ normal(1, 100);
  for(i in 1:Sobs){
    p[i] ~ beta(s1, s2) T[0.5/Nsites,];
    y[i] ~ binomial(Nsites, p[i]);
  }
}
generated quantities {
  real Sest;
  Sest = bruteforce_rng(Sobs, Smax, Nsites, Sobs, s1, s2);
}
