// // GOAL //
// // Model ONE lake at ONE grid cell with ONE time series //
// // 3.19.2026 //
// 
// // MAJOR MEMORY ISSUE (maybe resolved) //
// 
// // Imported from STEPPS - PREDICTION (2/6/2026)
// //
// // // Spatio-temporal vegetation model; veg maps predicted from pollen counts
// // // Latent vegetation modeled using a predictive process
// // // Linked to multinomial pollen counts through an additive log-ratio sum to one constraint
// //
// // weighting to take gamma = mean(gamma[k]) and w[i, j] instead of accounting for taxa
// 
data {

  // Data initiation
  // Only momdels k-1 taxa explicetly
  int<lower=0> K;       // number of species
  int<lower=0> N;       // number of cells
  int<lower=0> T;       // number of times
  int<lower=0> N_knots; // number of knots
  int<lower=0> N_cores; // number of cores

  array[N_cores * T, K] int y; // observed pollen counts per core-time

  // Taxon specific parameters
  // Choses other as a reference
  vector<lower=0>[K-1] rho;     // spatial covariance par (idk who it is k-1, it is taxa specific) HO-cpp starts counting vector length at 0 rather than 1
  vector<lower=0>[K-1] eta;     // space-time variance par

  // spatial mixing parameter (between local and regional pollen contribution)
  real<lower=0, upper=1> gamma; // local/long-distance dispersal par
 // vector<lower=0, upper=1>[K] gamma; // change to become taxa specific

  // real<lower=0> psi;            // pollen dispersal par - not actually used

  // Dirichelt precision parameter
  vector<lower=0>[K] phi;       // controls over dispersion

  array[N_cores] int idx_cores; // core cell indices

  // Exponential distance matricies
  matrix[N,N] d;                   // cell distance matrix
  matrix[N_knots,N_knots] d_knots; // knot distance matrix
  matrix[N,N_knots] d_inter;       // cells-knots distance matrix

  // pollen weighting matrix to see how much each cell contributes to other cell
  matrix[N_cores,N] w;             // weight matrix

  // NEVER USED
  matrix[T,T] lag;                 // time lag matrix

  int<lower=0> N_p;                // size of P in data file; will be 1 or N
  real P;                          // FIXME: want this to be a matrix OR a real
  
  //// ADD STAN ORIGINAL OUTPUTS TO IFNORM THE G LATENT FIELD
  int<lower=1> N_stepps;              // number of observations
  array[N_stepps] int cell_idx;      // which grid cell
  array[N_stepps] int time_idx;      // which time
  array[N_stepps] int taxa_idx;      // which taxa
  vector[N_stepps] stepps_mean;      // reference_prop
  vector<lower=0>[N_stepps] stepps_sd; // reference_sd
  real tau;
}

transformed data {
  int W;
  vector[K-1] eta2; //will be the space-time var parameter squared
  vector[N_knots] zeros; //Will be a vector of zeros

  matrix[N_knots,N_knots] Eye_knots; //Will be an identity matrix of dimensions N_knots X N_knots

  // Spatial covariance objects per taxon
  array[K-1] matrix[N_knots,N_knots] C_s;      // spatial covariance mat
  array[K-1] matrix[N_knots,N_knots] C_s_L;    // chol decomposition of C_s
  array[K-1] matrix[N_knots,N_knots] C_s_inv;  // inverse spatial covariance mat
  array[K-1] matrix[N,N_knots] c_s;            // spatial covariance mat

 //  matrix[N*T, N*T] M; // not used

  W = K-1;
  for (k in 1:W) eta2[k] = eta[k] * eta[k];
  for (i in 1:N_knots) zeros[i] = 0.0;

  //makes Eye_knots an identity matrix.
  for (i in 1:N_knots)
    for (j in 1:N_knots)
      Eye_knots[i,j] = (i==j ? 1.0 : 0.0);

  // FS - removed since never used and too large
  // for (j in 1:N*T)
  //   for (i in 1:N*T)
  //     M[i,j] = (i==j ? 1.0-P : -P);

  // construct spatial covariance matrix
  for (k in 1:W){
    // Construct covariance at knots
    C_s[k] = exp(-d_knots / rho[k]);

    // Add small jitter to diagonal for numerical stability (getting non-positive definite decomp)
    for (i in 1:N_knots)
      C_s[k][i,i] = C_s[k][i,i] + 1e-5;

    // Cholesky decomposition
    C_s_L[k] = cholesky_decompose(C_s[k]);

    // Inverse of the covariance
    C_s_inv[k]  = mdivide_right_tri_low(mdivide_left_tri_low(C_s_L[k], Eye_knots)', C_s_L[k])';

    // Cell-to-knot covariance
    c_s[k] = exp(-d_inter / rho[k]);
  }
}

parameters {

  // random walk smoothness parameter (0-20)
  real<lower=0, upper=20> ksi;

  vector<lower=0,upper=100>[W] sigma; // temporal varience HO - I am pretty sure this is SD
  vector<lower=0, upper=2>[W] lambda; // spactial range for temporal varience

  vector[W] mu;                        // scalar vector of species-varying overall adjustment term (mu_k in paper)
  array[W] vector[T] mu_t;             // time-varying mean - mu_{t,k}^t in Dawson ea 2019
  array[W] vector[N_knots] alpha_s;    // spatially-varying mean
  array[W*(T-1)] vector[N_knots] alpha_t; // temporal innovations
  array[W] vector[N*T] g;              // latent process (largest memory issue)

}

model {
  // declarations (storages)
  array[W] vector[N*T] mu_g;
  vector[N*T] sum_exp_g;
  array[N*T] vector[K] r; // Vegetation proportions

  // Covariance objects
  array[W] matrix[N, N_knots] q_s;
  array[W] matrix[N_knots, N_knots] Q_s;
  array[W] matrix[N_knots, N_knots] Q_s_L;
  array[W] matrix[N_knots, N_knots] Q_s_inv;

  vector[W] sigma2; //Gaussian proccess variance.

  real cvar;
  real qvar;

  vector[N] Halpha_s;
  vector[N] Halpha_t;
  array[T-1] vector[N] qQinv_alpha; // Probably nu_{., k}^2 (spatial process term) for knots

  row_vector[N_knots] c_i;//Variance for a given taxa and time
  row_vector[N_knots] q_i;
  vector[N*T] sqrtvar;

  // priors
  mu ~ normal(0,20);//Setting prior for species-specific overall adjustment term

  for (k in 1:W)
    mu_t[k][1] ~ normal(0.0, 20.0);//Set initial values for autoregressive mean ("time varying term")

  for (k in 1:W)//Set initial values for guassian process variance?
    sigma2[k] = sigma[k] * sigma[k];

  // innovations covariance
  for (k in 1:W){
    q_s[k]   = exp(-d_inter/lambda[k]);
    Q_s[k]   = exp(-d_knots/lambda[k]) + 1e-5; // add jitter
    Q_s_L[k] = cholesky_decompose(Q_s[k]);
    Q_s_inv[k] = mdivide_right_tri_low(mdivide_left_tri_low(Q_s_L[k], Eye_knots)', Q_s_L[k])';
  }

  // LATENT GP (gaussian process)

  for (k in 1:W){
    // spatially-varying mean
    alpha_s[k] ~ multi_normal_cholesky(zeros, eta2[k] * C_s_L[k]);

    // Predictive process projection to cells
    Halpha_s = c_s[k] * (C_s_inv[k] * alpha_s[k]);

    for (i in 1:N)
      for (t in 1:T)
        mu_g[k][(i-1)*T + t] = Halpha_s[i];

    // orthogonal decomposition
    mu_g[k] = mu_g[k];

   //  print("Here");
   // print("k=", k);

    // temporal innovations - (random walk)
    alpha_t[(k-1)*(T-1) + 1] ~ multi_normal_cholesky(zeros, cholesky_decompose(1/sigma2[k] * Q_s_inv[k]));

    for (t in 2:(T-1))
      alpha_t[(k-1)*(T-1) + t] ~ multi_normal_cholesky(alpha_t[(k-1)*(T-1)+t-1], cholesky_decompose(1/sigma2[k] * Q_s_inv[k]));

    for (t in 1:(T-1))
      qQinv_alpha[t] = q_s[k] * Q_s_inv[k] * alpha_t[(k-1)*(T-1)+t]; // Spatial process term for all time points.

    // time-varying mean (random walk)
    for (i in 2:T)
      mu_t[k][i] ~ normal(mu_t[k][i-1], ksi);

    // builds full spatio-temporal mean and variance
    for (i in 1:N){
      c_i = row(c_s[k], i);//variance between knots and cells for a given taxa and time
      cvar = eta2[k] * c_i * C_s_inv[k] * c_i'; //Maybe the covariance for the mu terms in the gaussian process? Not sure.

      q_i = row(q_s[k], i); //Variance of innovations for a given taxa and cell?
      qvar = sigma2[k] * q_i * Q_s_inv[k] * q_i';

      // t = 1
      mu_g[k][(i-1)*T + 1] = mu[k] + mu_t[k][1] + mu_g[k][(i-1)*T + 1];
      //sqrtvar[(i-1)*T + 1] = sqrt(eta2[k] - cvar); # might be causing converging issues
      sqrtvar[(i-1)*T + 1] = sqrt(fmax(eta2[k] - cvar, 1e-10)); // take max

      // t > 1
      for (t in 2:T){
        mu_g[k][(i-1)*T + t] = mu[k] + mu_t[k][t] + mu_g[k][(i-1)*T + t] + qQinv_alpha[t-1][i];
        sqrtvar[(i-1)*T + t] = sqrt(eta2[k] - cvar + sigma2[k] - qvar);
      }
    }

    // Sample latent Gaussian field
    for (i in 1:N*T)
      g[k][i] ~ normal(mu_g[k][i], sqrtvar[i]);
  }

  // sum exponential of process vals
  for (i in 1:N*T){
    sum_exp_g[i] = 0.0;
    for (k in 1:W)
      sum_exp_g[i] += exp(g[k,i]);
  }

  for (k in 1:W)
    for (i in 1:N*T)
      r[i,k] = exp(g[k,i]) / (1 + sum_exp_g[i]);

  for (i in 1:N*T)
    r[i,K] = 1 / (1 + sum_exp_g[i]);
  
  ///// FINLEY - ADD THE STEPPS ORIGNAL OUTPUT ////////////  
  for (n in 1:N_stepps) {
  int i = cell_idx[n];
  int t = time_idx[n];
  int k = taxa_idx[n];

  int idx = (i-1)*T + t;
  
  ////// tau = tuning parameter
  stepps_mean[n] ~ normal(r[idx,k], tau * stepps_sd[n]);
  }
  ///// FINLEY - END ADD THE STEPPS ORIGNAL OUTPUT //////////// 
  
  // link to pollen (pollen dispersal + likelihood)
  {
    array[N_cores*T] vector[K] r_new;
    vector[K] out_sum;
    real sum_w;
    int idx_core;
    real N_grains;
    real A;
    vector[K] kappa;

    for (i in 1:N_cores){
      for (t in 1:T){
        idx_core = (idx_cores[i]-1)*T + t;

        // local vegetation contribution
        r_new[(i-1)*T + t] = gamma * r[idx_core]; // old has gamma = scalar
        // r_new[(i-1)*T + t] = gamma[K] * r[idx_core];

        // long-distance pollen contribution
        out_sum = rep_vector(0.0, K);
        sum_w = 0;

        for (j in 1:N)
          if (d[idx_cores[i],j] > 0)
            out_sum += w[i,j] * r[(j-1)*T + t];

        sum_w = sum(out_sum);

        // end long-distance pollen contribution

        // if there is long distance pollen contribution:
        if (sum_w > 0)
          r_new[(i-1)*T + t] += (1-gamma) * out_sum / sum_w; // old has gamma = scalar
         // r_new[(i-1)*T + t] += (1-gamma[K]) * out_sum / sum_w;
      }
    }

    // Dirichlet–Multinomial likelihood
    for (i in 1:N_cores*T){
      if (sum(y[i]) > 0){
        kappa = phi .* r_new[i];
        A = sum(kappa);
        N_grains = sum(y[i]);

        target += lgamma(N_grains + 1) + lgamma(A) - lgamma(N_grains + A);

        for (k in 1:K)
          target += -lgamma(y[i,k] + 1) + lgamma(y[i,k] + kappa[k]) - lgamma(kappa[k]);
      }
    }
  }
}

// save additional outputs including observed pollen
generated quantities {
  array[N*T] vector[K] r;
  array[N_cores*T] vector[K] r_new;

  vector[N*T] sum_exp_g;

  // rebuild r exactly as in model
  for (i in 1:N*T) {
    sum_exp_g[i] = 0;

    for (k in 1:W)
      sum_exp_g[i] += exp(g[k,i]);

    for (k in 1:W)
      r[i,k] = exp(g[k,i]) / (1 + sum_exp_g[i]);

    r[i,K] = 1 / (1 + sum_exp_g[i]);
  }

  //  compute r_new safely
  for (i in 1:N_cores) {
    for (t in 1:T) {

      int idx_core = (idx_cores[i]-1)*T + t;

      r_new[(i-1)*T + t] = gamma * r[idx_core];

      vector[K] out_sum = rep_vector(0.0, K);
      real sum_w = 0;

      for (j in 1:N)
        if (d[idx_cores[i], j] > 0)
          out_sum += w[i,j] * r[(j-1)*T + t];

      sum_w = sum(out_sum);

      if (sum_w > 0)
        r_new[(i-1)*T + t] += (1 - gamma) * out_sum / sum_w;
    }
  }
}



// // GOAL:
// // Single lake, single cell, time series only
// // Spatial structure REMOVED
// 
// data {
//   int<lower=1> K;     // number of taxa (12)
//   int<lower=1> T;     // number of time points (20)
//   array[T, K] int y;  // pollen counts per time
//   vector<lower=0>[K] phi;  // Dirichlet precision
// }
// 
// transformed data {
//   int W = K - 1;  // controls over dispersion (but there won't be dispersion...?)
// }
// 
// parameters {
//   vector[W] mu;
// 
//   vector<lower=0>[W] sigma_rw;
//   vector<lower=0>[W] sigma_g;
// 
//   array[W] vector[T] mu_t;
//   array[W] vector[T] g;
// }
// 
// 
// // actual model
// model {
//   vector[T] sum_exp_g;
//   array[T] vector[K] r;
// 
//   mu ~ normal(0, 20);
//   sigma_rw ~ normal(0, 2);
//   sigma_g ~ normal(0, 2);
// 
//   for (k in 1:W) {
//     mu_t[k][1] ~ normal(0, 5);
// 
//     for (t in 2:T)
//       mu_t[k][t] ~ normal(mu_t[k][t-1], sigma_rw[k]);
// 
//     for (t in 1:T)
//       g[k][t] ~ normal(mu[k] + mu_t[k][t], sigma_g[k]);
//   }
// 
//   //  LATENT PROCESS
//   for (k in 1:W) {
//     for (t in 1:T) {
//       // g[k][t] ~ normal(mu[k] + mu_t[k][t], sigma[k]);
//       g[k][t] ~ normal(mu[k] + mu_t[k][t], sigma[k]);
//     }
//   }
// 
//   //
//   for (t in 1:T) {
//     sum_exp_g[t] = 0;
// 
//     for (k in 1:W)
//       sum_exp_g[t] += exp(g[k][t]);
// 
//     for (k in 1:W)
//       r[t][k] = exp(g[k][t]) / (1 + sum_exp_g[t]);
// 
//     r[t][K] = 1 / (1 + sum_exp_g[t]);
//   }
// 
//   // DIRICHLET-MULTINOMIAL 
//   for (t in 1:T) {
//     if (sum(y[t]) > 0) {
//       vector[K] kappa = phi .* r[t];
//       real A = sum(kappa);
//       real N = sum(y[t]);
// 
//       target += lgamma(N + 1) + lgamma(A) - lgamma(N + A);
// 
//       for (k in 1:K)
//         target += -lgamma(y[t,k] + 1)
//                   + lgamma(y[t,k] + kappa[k])
//                   - lgamma(kappa[k]);
//     }
//   }
// }

