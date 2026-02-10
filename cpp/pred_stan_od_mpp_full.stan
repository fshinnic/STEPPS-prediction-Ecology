// MAJOR MEMORY ISSUE //

// Imported from STEPPS - PREDICTION (2/6/2026)
//
// // Spatio-temporal vegetation model; veg maps predicted from pollen counts
// // Latent vegetation modeled using a predictive process
// // Linked to multinomial pollen counts through an additive log-ratio sum to one constraint
//
// weighting to take gamma = mean(gamma[k]) and w[i, j] instead of accounting for taxa

data {

  int<lower=0> K;       // number of species
  int<lower=0> N;       // number of cells
  int<lower=0> T;       // number of times
  int<lower=0> N_knots; // number of knots
  int<lower=0> N_cores; // number of cores

  array[N_cores * T, K] int y;

  vector<lower=0>[K-1] rho;     // spatial covariance par (idk who it is k-1, it is taxa specific)
  vector<lower=0>[K-1] eta;     // space-time variance par


   real<lower=0, upper=1> gamma; // local/long-distance dispersal par
 // vector<lower=0, upper=1>[K] gamma; // change to become taxa specific

 // real<lower=0> psi;            // pollen dispersal par - not actually used

  vector<lower=0>[K] phi;       // dirichlet precision par

  array[N_cores] int idx_cores; // core cell indices

  matrix[N,N] d;                   // cell distance matrix
  matrix[N_knots,N_knots] d_knots; // knot distance matrix
  matrix[N,N_knots] d_inter;       // cells-knots distance matrix

  matrix[N_cores,N] w;             // weight matrix

  matrix[T,T] lag;                 // time lag matrix

  int<lower=0> N_p;                // size of P in data file; will be 1 or N
  real P;                          // FIXME: want this to be a matrix OR a real
}

transformed data {
  int W;
  vector[K-1] eta2;
  vector[N_knots] zeros;

  matrix[N_knots,N_knots] Eye_knots;

  // modern array syntax
  array[K-1] matrix[N_knots,N_knots] C_s;      // spatial covariance mat
  array[K-1] matrix[N_knots,N_knots] C_s_L;    // chol decomposition of C_s
  array[K-1] matrix[N_knots,N_knots] C_s_inv;  // inverse spatial covariance mat
  array[K-1] matrix[N,N_knots] c_s;            // spatial covariance mat

  matrix[N*T, N*T] M;

  W = K-1;
  for (k in 1:W) eta2[k] = eta[k] * eta[k];
  for (i in 1:N_knots) zeros[i] = 0.0;

  for (i in 1:N_knots)
    for (j in 1:N_knots)
      Eye_knots[i,j] = (i==j ? 1.0 : 0.0);

  // FS - removed since never used and too large
  // for (j in 1:N*T)
  //   for (i in 1:N*T)
  //     M[i,j] = (i==j ? 1.0-P : -P);

  // construct spatial covariance matrix
  for (k in 1:W){
    C_s[k]   = exp(-d_knots/rho[k]);
    C_s_L[k] = cholesky_decompose(C_s[k]);
    C_s_inv[k]  = mdivide_right_tri_low(mdivide_left_tri_low(C_s_L[k], Eye_knots)', C_s_L[k])';
    c_s[k]      = exp(-d_inter/rho[k]);
  }
}

parameters {

  real<lower=0, upper=20> ksi;

  vector<lower=0,upper=100>[W] sigma;
  vector<lower=0, upper=2>[W] lambda;

  vector[W] mu;                        // scalar vector
  array[W] vector[T] mu_t;             // time-varying mean
  array[W] vector[N_knots] alpha_s;    // spatially-varying mean
  array[W*(T-1)] vector[N_knots] alpha_t; // temporal innovations
  array[W] vector[N*T] g;              // latent process

}

model {
  // declarations
  array[W] vector[N*T] mu_g;
  vector[N*T] sum_exp_g;
  array[N*T] vector[K] r;

  array[W] matrix[N, N_knots] q_s;
  array[W] matrix[N_knots, N_knots] Q_s;
  array[W] matrix[N_knots, N_knots] Q_s_L;
  array[W] matrix[N_knots, N_knots] Q_s_inv;

  vector[W] sigma2;

  real cvar;
  real qvar;

  vector[N] Halpha_s;
  vector[N] Halpha_t;
  array[T-1] vector[N] qQinv_alpha;

  row_vector[N_knots] c_i;
  row_vector[N_knots] q_i;
  vector[N*T] sqrtvar;

  // priors
  mu ~ normal(0,20);

  for (k in 1:W)
    mu_t[k][1] ~ normal(0.0, 20.0);

  for (k in 1:W)
    sigma2[k] = sigma[k] * sigma[k];

  // innovations covariance
  for (k in 1:W){
    q_s[k]   = exp(-d_inter/lambda[k]);
    Q_s[k]   = exp(-d_knots/lambda[k]);
    Q_s_L[k] = cholesky_decompose(Q_s[k]);
    Q_s_inv[k] = mdivide_right_tri_low(mdivide_left_tri_low(Q_s_L[k], Eye_knots)', Q_s_L[k])';
  }

  for (k in 1:W){
    // spatially-varying mean
    alpha_s[k] ~ multi_normal_cholesky(zeros, eta2[k] * C_s_L[k]);

    Halpha_s = c_s[k] * (C_s_inv[k] * alpha_s[k]);
    for (i in 1:N)
      for (t in 1:T)
        mu_g[k][(i-1)*T + t] = Halpha_s[i];

    // orthogonal decomposition
    mu_g[k] = mu_g[k];

   //  print("Here");
   // print("k=", k);

    // temporal innovations
    alpha_t[(k-1)*(T-1) + 1] ~ multi_normal_cholesky(zeros, cholesky_decompose(1/sigma2[k] * Q_s_inv[k]));
    for (t in 2:(T-1))
      alpha_t[(k-1)*(T-1) + t] ~ multi_normal_cholesky(alpha_t[(k-1)*(T-1)+t-1], cholesky_decompose(1/sigma2[k] * Q_s_inv[k]));

    for (t in 1:(T-1))
      qQinv_alpha[t] = q_s[k] * Q_s_inv[k] * alpha_t[(k-1)*(T-1)+t];

    // time-varying mean
    for (i in 2:T)
      mu_t[k][i] ~ normal(mu_t[k][i-1], ksi);

    for (i in 1:N){
      c_i = row(c_s[k], i);
      cvar = eta2[k] * c_i * C_s_inv[k] * c_i';

      q_i = row(q_s[k], i);
      qvar = sigma2[k] * q_i * Q_s_inv[k] * q_i';

      mu_g[k][(i-1)*T + 1] = mu[k] + mu_t[k][1] + mu_g[k][(i-1)*T + 1];
      sqrtvar[(i-1)*T + 1] = sqrt(eta2[k] - cvar);

      for (t in 2:T){
        mu_g[k][(i-1)*T + t] = mu[k] + mu_t[k][t] + mu_g[k][(i-1)*T + t] + qQinv_alpha[t-1][i];
        sqrtvar[(i-1)*T + t] = sqrt(eta2[k] - cvar + sigma2[k] - qvar);
      }
    }

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

  // link to pollen
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

        r_new[(i-1)*T + t] = gamma * r[idx_core]; // old has gamma = scalar
        // r_new[(i-1)*T + t] = gamma[K] * r[idx_core];

        out_sum = rep_vector(0.0, K);
        sum_w = 0;

        for (j in 1:N)
          if (d[idx_cores[i],j] > 0)
            out_sum += w[i,j] * r[(j-1)*T + t];

        sum_w = sum(out_sum);

        if (sum_w > 0)
          r_new[(i-1)*T + t] += (1-gamma) * out_sum / sum_w; // old has gamma = scalar
         // r_new[(i-1)*T + t] += (1-gamma[K]) * out_sum / sum_w;
      }
    }

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




// =========================
// MEMORY ALLOCATION AUDIT
// =========================
// 
// // ---------- data block ----------
// data {
//   int<lower=0> K;                 // small
//   int<lower=0> N;                 // large driver
//   int<lower=0> T;                 // moderate
//   int<lower=0> N_knots;           // moderate
//   int<lower=0> N_cores;           // moderate
// 
//   array[N_cores * T, K] int y;    // O(N_cores * T * K)
// 
//   vector[K-1] rho;                // small
//   vector[K-1] eta;                // small
//   real gamma;                     // scalar
//   vector[K] phi;                  // small
// 
//   array[N_cores] int idx_cores;   // small
// 
//   matrix[N,N] d;                  //  O(N^2) VERY LARGE
//   matrix[N_knots,N_knots] d_knots;// O(N_knots^2)
//   matrix[N,N_knots] d_inter;      // O(N * N_knots)
// 
//   matrix[N_cores,N] w;            // O(N_cores * N)
// 
//   matrix[T,T] lag;                // small (unused)
//   int N_p;                        // scalar
//   real P;                         // scalar
// }
// 
// 
// // ---------- transformed data ----------
// transformed data {
//   int W;                          // scalar
//   vector[K-1] eta2;               // small
//   vector[N_knots] zeros;          // O(N_knots)
// 
//   matrix[N_knots,N_knots] Eye_knots;  // O(N_knots^2)
// 
//   // Per-taxon spatial objects
//   array[K-1] matrix[N_knots,N_knots] C_s;      // O(W * N_knots^2)
//   array[K-1] matrix[N_knots,N_knots] C_s_L;    // O(W * N_knots^2)
//   array[K-1] matrix[N_knots,N_knots] C_s_inv;  //  O(W * N_knots^2) + autodiff
//   array[K-1] matrix[N,N_knots] c_s;             //  O(W * N * N_knots)
// 
//   matrix[N*T, N*T] M;             // O((N*T)^2) — EVEN IF UNUSED
// }
// 
// 
// // ---------- parameters ----------
// parameters {
//   real ksi;                       // scalar
// 
//   vector[W] sigma;                // small
//   vector[W] lambda;               // small
//   vector[W] mu;                   // small
// 
//   array[W] vector[T] mu_t;        // O(W * T)
//   array[W] vector[N_knots] alpha_s; // O(W * N_knots)
//   array[W*(T-1)] vector[N_knots] alpha_t; // O(W * T * N_knots)
// 
//   array[W] vector[N*T] g;         //  O(W * N * T) HUGE + autodiff
// }
// 
// 
// // ---------- model (local allocations) ----------
// model {
//   // Latent mean storage
//   array[W] vector[N*T] mu_g;      //  O(W * N * T)
//   vector[N*T] sum_exp_g;          // O(N * T)
// 
//   // Probability simplex for *all cells*
//   array[N*T] vector[K] r;         //  O(N * T * K)
// 
//   // Predictive-process covariance
//   array[W] matrix[N, N_knots] q_s;       //  O(W * N * N_knots)
//   array[W] matrix[N_knots,N_knots] Q_s;  // O(W * N_knots^2)
//   array[W] matrix[N_knots,N_knots] Q_s_L;// O(W * N_knots^2)
//   array[W] matrix[N_knots,N_knots] Q_s_inv;//  heavy autodiff
// 
//   vector[W] sigma2;               // small
// 
//   // Temporary per-cell buffers
//   vector[N] Halpha_s;             // O(N)
//   vector[N] Halpha_t;             // O(N)
//   array[T-1] vector[N] qQinv_alpha; //  O(N * T)
// 
//   row_vector[N_knots] c_i;        // O(N_knots)
//   row_vector[N_knots] q_i;        // O(N_knots)
// 
//   vector[N*T] sqrtvar;            // O(N * T)
// 
//   // Likelihood-only storage (reasonable)
//   array[N_cores*T] vector[K] r_new; // O(N_cores * T * K)
// }
// 
