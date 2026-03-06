data {
  int<lower=1> N_train;                         // Number of training observations
  int<lower=1> N_test;                          // Number of test observations
  int<lower=1> R;                               // Number of spatial units (regions)
  int<lower=1> T;                               // Number of time points (days)

  array[N_train] int<lower=1, upper=R> region_train;   // Region index for training observations
  array[N_train] int<lower=1, upper=T> time_train;     // Time index for training observations
  array[N_train] int<lower=1,upper=7> dow_train;       // Day-of-week index (1–7)
  array[N_train] real y_log_train;                     // Log-transformed response variable

  array[N_test] int<lower=1, upper=R> region_test;     // Region index for test observations
  array[N_test] int<lower=1, upper=T> time_test;       // Time index for test observations
  array[N_test] int<lower=1,upper=7> dow_test;         // Day-of-week index for test set

  vector[R] pop;                                // Population size for each region
  matrix[R, R] Dmat_space;                      // Pairwise geographic distance matrix
}

transformed data {
  real phi_s_min;
  real phi_s_max;
  matrix[R, R] Dmat_space_no_diag;

  // Lower bound for spatial decay parameter
  phi_s_min = 3 / max(Dmat_space);

  // Replace diagonal entries to avoid zero distances
  Dmat_space_no_diag = Dmat_space;
  for (i in 1:R)
    Dmat_space_no_diag[i, i] = positive_infinity();

  // Upper bound for spatial decay parameter
  phi_s_max = 3 / min(to_array_1d(Dmat_space_no_diag));
}

parameters {
  real mu;                                      // Global intercept
  real<lower=0> sigma2;                         // Variance of the spatiotemporal latent process

  real<lower=-0.99, upper=0.99> alpha;          // AR(1) temporal autocorrelation parameter
  matrix[R, T] w_t;                             // Temporal latent effects (AR(1) process per region)

  vector[R] w_s;                                // Spatial latent effects

  real<lower=phi_s_min, upper = phi_s_max> phi_s;  // Spatial decay parameter
  real<lower=0, upper=1> rho;                      // Weight balancing spatial vs temporal effects

  vector[6] gamma;                               // Day-of-week fixed effects (centered)
  real<lower=0> sigma2_y;                         // Observation noise variance
  real beta;                                      // Population effect coefficient
}

transformed parameters {

  // Extend day-of-week effects with baseline constraint
  vector[7] gamma_ext;
  for (d in 1:6)
    gamma_ext[d] = gamma[d];
  gamma_ext[7] = 0;

  real<lower = 0> sigma_y = sqrt(sigma2_y);     // Observation noise standard deviation
  real<lower=0> sigma_t = sqrt(1-alpha^2);      // Innovation std for AR(1) process
  real<lower=0> sigma = sqrt(sigma2);           // Std of latent spatiotemporal process

  // Spatial covariance matrix (exponential kernel)
  matrix[R, R] Sigma_cov;
  for (i in 1:R)
    for (j in 1:R)
      Sigma_cov[i, j] = exp(-phi_s * Dmat_space[i, j]);
}

model {

  // Priors
  mu ~ normal(0, 100);
  sigma2 ~ inv_gamma(1, 1);

  alpha ~ uniform(-1, 1);
  phi_s ~ uniform(phi_s_min, phi_s_max);

  rho ~ beta(1, 2);

  gamma ~ normal(0, 10);

  sigma2_y ~ inv_gamma(1, 1);

  beta ~ normal(0, 1000);

  // Spatial Gaussian process
  w_s ~ multi_normal(rep_vector(0, R), Sigma_cov);

  // Temporal AR(1) process for each region
  for (r in 1:R) {

    // Initial state drawn from stationary distribution
    w_t[r, 1] ~ normal(0, sigma_t / fmax(sqrt(1 - square(alpha)), 1e-6));

    for (t in 2:T)
      w_t[r, t] ~ normal(alpha * w_t[r, t - 1], sigma_t);
  }

  // Likelihood
  for (n in 1:N_train) {

    int s = region_train[n];
    int t = time_train[n];
    int d = dow_train[n];

    real mu_y =
      mu
      + gamma_ext[d]
      + sigma * (rho * w_s[s] + (1 - rho) * w_t[s, t])
      + beta * log(pop[s]);

    y_log_train[n] ~ normal(mu_y, sigma_y);
  }
}

generated quantities {

  array[N_test] real y_log_pred;
  vector[N_train] log_lik;

  for (n in 1:N_test) {

    int s = region_test[n];
    int t = time_test[n];
    int d = dow_test[n];

    real mu_y =
      mu
      + gamma_ext[d]
      + sigma * (rho * w_s[s] + (1 - rho) * w_t[s, t])
      + beta * log(pop[s]);

    y_log_pred[n] = normal_rng(mu_y, sigma_y);
  }

  for (n in 1:N_train) {

    int s = region_train[n];
    int t = time_train[n];
    int d = dow_train[n];

    real mu_y =
      mu
      + gamma_ext[d]
      + sigma * (rho * w_s[s] + (1 - rho) * w_t[s, t])
      + beta * log(pop[s]);

    log_lik[n] = normal_lpdf(y_log_train[n] | mu_y, sigma_y);
  }
}
