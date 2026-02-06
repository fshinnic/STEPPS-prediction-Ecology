// Improted from STEPPS - PREDICTION (2/6/2026)

// Spatio-temporal vegetation model; veg maps predicted from pollen counts
// Latent vegetation modelled using a predictive process
// Linked to multinomial pollen counts through an additive log-ratio sum to one constraint

data {

  int<lower=0> K;       // number of species
  int<lower=0> N;       // number of cells
  int<lower=0> T;       // number of times
  int<lower=0> N_knots; // number of knots
  int<lower=0> N_cores; // number of cores 

  int y[N_cores*T,K];
  
  vector<lower=0>[K-1] rho;     // spatial covariance par
  vector<lower=0>[K-1] eta;     // space-time variance par
  real<lower=0, upper=1> gamma; // local/long-distance dispersal par 
  real<lower=0> psi;            // pollen dispersal par

  vector<lower=0>[K] phi;     // dirichlet precision par
  
  int idx_cores[N_cores];          // core cell indices

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
  
  matrix[N_knots,N_knots]  Eye_knots;
  
  matrix[N_knots,N_knots] C_s[K-1]; // spatial covariance mat
  matrix[N_knots,N_knots] C_s_L[K-1]; // chol decomposition of C_s
  matrix[N_knots,N_knots] C_s_inv[K-1]; // inverse spatial covariance mat
  matrix[N,N_knots] c_s[K-1]; // spatial covariance mat
  
  matrix[N*T, N*T] M;
  
  W      <- K-1;
  for (k in 1:W) eta2[k]   <- eta[k] * eta[k];
  for (i in 1:N_knots) zeros[i] <- 0.0;

  for (i in 1:N_knots)
    for (j in 1:N_knots)
      if (i == j) Eye_knots[i,j] <- 1.0;
      else        Eye_knots[i,j] <- 0;

  for (j in 1:N*T) 
    for (i in 1:N*T) 
      if (i==j) M[i,j] <- 1.0-P;
      else      M[i,j] <- -P;
      
  // construct spatial covariance matrix
  for (k in 1:W){
    C_s[k]   <- exp(-d_knots/rho[k]);
    C_s_L[k] <- cholesky_decompose(C_s[k]);
    C_s_inv[k]  <- mdivide_right_tri_low(mdivide_left_tri_low(C_s_L[k], Eye_knots)', C_s_L[k])'; 
    c_s[k]      <- exp(-d_inter/rho[k]);
  } 
}

parameters {

  real<lower=0, upper=20> ksi;   

  vector<lower=0,upper=100>[W] sigma;
  vector<lower=0, upper=2>[W] lambda;

  vector[W] mu;
  vector[T] mu_t[W];
  vector[N_knots] alpha_s[W];
  vector[N_knots] alpha_t[W*(T-1)];
  vector[N*T] g[W];
}

model {

  //declarations
  vector[N*T] mu_g[W];
  vector[N*T] sum_exp_g;  
  vector[K]   r[N*T];     // cell proportions
  
  matrix[N, N_knots] q_s[W];
  matrix[N_knots, N_knots] Q_s[W];
  matrix[N_knots, N_knots] Q_s_L[W];
  matrix[N_knots, N_knots] Q_s_inv[W];

  vector[W] sigma2;

  real cvar; 
  real qvar;

  vector[N] Halpha_s; 
  vector[N] Halpha_t;
  vector[N] qQinv_alpha[T-1];
 
  row_vector[N_knots] c_i; 
  row_vector[N_knots] q_i; 
  vector[N*T] sqrtvar;

  // priors 
  mu       ~ normal(0,20);
  
  for (k in 1:W){
    //    mu_t[k][1] ~ normal(mu[k], sqrt((ksi * ksi)/(1 - omega * omega)) ); 
    mu_t[k][1] ~ normal(0.0, 20.0); 
  }

  for (k in 1:W){
    sigma2[k] <- sigma[k] * sigma[k];  
  }
 
  // innovations covariance
  for (k in 1:W){
    q_s[k]      <- exp(-d_inter/lambda[k]);
    Q_s[k]      <- exp(-d_knots/lambda[k]);
    Q_s_L[k]    <- cholesky_decompose(Q_s[k]);
    Q_s_inv[k]  <- mdivide_right_tri_low(mdivide_left_tri_low(Q_s_L[k], Eye_knots)', Q_s_L[k])'; 
  }

  for (k in 1:W){ 

    //    print("sqrtvar : ", dims(sqrtvar));


    // spatially-varying mean
    // alpha_s[k] ~ multi_normal_prec(zeros, 1/eta2[k] * C_s_inv[k]);  
    alpha_s[k] ~ multi_normal_cholesky(zeros, eta2[k] * C_s_L[k]);  
    
    Halpha_s <- c_s[k] * (C_s_inv[k] * alpha_s[k]);
    for (i in 1:N){
      for (t in 1:T){
        mu_g[k][(i-1) * T + t] <- Halpha_s[i];
      }
    } 

    // orthogonal decomposition
    //mu_g[k]  <- M * mu_g[k];
    mu_g[k]  <- mu_g[k];

    print("Here")

    // innovations
    // really the second, but the first are set to zeros
    // having issues parsing parenthesses
    alpha_t[((k-1)*(T-1) + 1)] ~ multi_normal_prec(zeros, 1/sigma2[k] * Q_s_inv[k]);


    // get string parsing error
    
    print("k = " + to_string(k));


    for (t in 2:(T-1)){
      print("index = ", (k-1)*(T-1)+t-1);
      alpha_t[(k-1)*(T-1)+t] ~ multi_normal_prec(alpha_t[(k-1)*(T-1)+t-1], 1/sigma2[k] *  Q_s_inv[k]); 
    }

    for (t in 1:(T-1)){
      qQinv_alpha[t] <- q_s[k] * Q_s_inv[k] * alpha_t[(k-1)*(T-1)+t];
    }

    // time-varying mean
    for (i in 2:T){
      mu_t[k][i] ~ normal(mu_t[k][i-1], ksi);
    }    

    for (i in 1:N){

      c_i  <- row(c_s[k], i);
      cvar <- eta2[k] * c_i * C_s_inv[k] * c_i';

      q_i  <- row(q_s[k], i);
      qvar <- sigma2[k] * q_i * Q_s_inv[k] * q_i';
      
      mu_g[k][(i-1) * T + 1]  <- mu[k] + mu_t[k][1] + mu_g[k][(i-1) * T + 1];
      sqrtvar[(i-1)*T + 1]    <- sqrt(eta2[k] - cvar);

      for (t in 2:T){
	mu_g[k][(i-1) * T + t]  <- mu[k] + mu_t[k][t] + mu_g[k][(i-1) * T + t] +  qQinv_alpha[t-1][i];
	sqrtvar[(i-1)*T + t]    <- sqrt(eta2[k] - cvar + sigma2[k] - qvar);	
      }
    }
    
    for (i in 1:N*T){
      g[k][i] ~ normal(mu_g[k][i], sqrtvar[i]);
    }

  }

  // sum exponential of process vals for W taxa by cell
  for (i in 1:N*T) {
    sum_exp_g[i] <- 0.0;
    for (k in 1:W)
      sum_exp_g[i] <- sum_exp_g[i] + exp(g[k,i]);
  }

  for (k in 1:W)
    for (i in 1:N*T)
      r[i,k] <- exp(g[k,i]) / (1 + sum_exp_g[i]);

  for (i in 1:N*T)
      r[i,K] <- 1 / (1 + sum_exp_g[i]);
      
  
  { // link to pollen! 
  vector[K] r_new[N_cores*T];
  vector[K] out_sum;    
  real sum_w;
  int idx_core;
  
  for (i in 1:N_cores){
    for (t in 1:T){    
      idx_core <- (idx_cores[i]-1)*T+t;
      r_new[(i-1)*T + t] <- gamma*r[idx_core];

      for (k in 1:K) {out_sum[k] <- 0;}
      sum_w <- 0; 

      for (j in 1:N){
        if (d[idx_cores[i],j] > 0){ 
          out_sum <- out_sum + w[i,j]*r[(j-1)*T+t];
        }  
       }
     
      sum_w   <- sum(out_sum);
    
      r_new[(i-1)*T + t] <- r_new[(i-1)*T + t] + out_sum*(1-gamma)/sum_w;   
  // link composition vector to count data through multinomial
    }
  }
  
  {
   real N_grains;
   real A;
   vector[K] kappa;

   for (i in 1:N_cores*T){ 
      
     if (sum(y[i]) > 0){ 
      kappa <- phi .* r_new[i];
      
      A <- sum(kappa);
      N_grains <- sum(y[i]);
      
      //lp__ <- lp__ + lgamma(N_grains + 1) + lgamma(A) - lgamma(N_grains + A);
      
      increment_log_prob(lgamma(N_grains + 1) + lgamma(A) - lgamma(N_grains + A));
      //for (k in 1:K) lp__ <- lp__ - lgamma(y[i,k] + 1) + lgamma(y[i,k] + kappa[k]) - lgamma(kappa[k]);
      
      for (k in 1:K) increment_log_prob(- lgamma(y[i,k] + 1) + lgamma(y[i,k] + kappa[k]) - lgamma(kappa[k]));
     }
    
    } 
  } // end scope
      

  } // end scope
}
