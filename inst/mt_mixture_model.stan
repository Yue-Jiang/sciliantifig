data {
  int<lower = 0> N;
  int<lower = 0> y[N];
  int n_groups;
}
parameters {
  //ordered doesn't work with constraints
  //so set constraints later in transformed
  ordered[n_groups] p_logit;
  simplex[n_groups] Theta;
}
transformed parameters {
  vector<lower=0, upper=1>[n_groups] p = inv_logit(p_logit);
}
model {
  vector[n_groups] contributions;
  // priors
  p ~ beta(5, 5);
  Theta ~ dirichlet(rep_vector(1.0, n_groups));

  // likelihood
  for(i in 1:N) {
    for(k in 1:n_groups) {
      contributions[k] = log(Theta[k]) + binomial_lpmf(y[i] | 19, p[k]);
    }
    target += log_sum_exp(contributions);
  }
}

// data {
//  int<lower = 0> N;
//  int<lower=0> y[N];
// }
//
// parameters {
//   simplex[3] theta;
//   real<lower=0, upper=1> p1;
//   real<lower=0, upper=1> p2;
//   real<lower=0, upper=1> p3;
// }
//
// model {
//  p1 ~ beta(0.9, 0.1);
//  p2 ~ beta(0.5, 0.5);
//  p3 ~ beta(0.1, 0.9);
//  for (n in 1:N)
//    target += log_mix(theta,
//                      binomial_lpmf(y[n] | 19, p1),
//                      binomial_lpmf(y[n] | 19, p2),
//                      binomial_lpmf(y[n] | 19, p3));
// }
