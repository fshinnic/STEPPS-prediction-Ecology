# instal package to run .stan from Rstudio
library(cmdstanr)
#cmdstanr::install_cmdstan()

# version check
cmdstanr::cmdstan_version()
# 2.38.0

############## required parameters #################
# K = taxa
# N =   
# T =  
# N_knots 
# N_cores
# y , matrix [N_cores*T, K]
#   rho                 # vector [K-1] > 0
#   eta            # vector [K-1] > 0
#   gamma          # real in [0,1]
#   psi         # real > 0 (declared, unused)
#   phi             # vector [K] > 0
#   idx_cores      # int vector [N_cores]
#   d                   # matrix [N, N]
#   d_knots      # matrix [N_knots, N_knots]
#   d_inter      # matrix [N, N_knots]
#   w                      # matrix [N_cores, N]
#   lag                # matrix [T, T] (declared, unused)
#   N_p _                 # int (declared, unused)
#   P                      # real scalar

# FS chosen files
fname <- file.path(dirName, subDir, 'input')  # same as used when saving
rdata_file <- paste0(fname, '.rdata')

# Load the file
load(rdata_file)

####### input parameters #######
# K, N, T, N_cores, N_knots, res, gamma, phi, rho, eta, y, idx_cores, 
# d_knots, d_inter, w, d, lag, P, N_p, meta_pol, meta_pol_all, knot_coords, 
# centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls

# missing the psi parameter, but psi not actually used in the pred_stan_od_mpp_full.stan file

###### run stan? #####

# Compile the model
stan_file <- file.path("cpp/pred_stan_od_mpp_full.stan")
mod <- cmdstan_model(stan_file) # now works!

# Run the model
fit <- mod$sample(
  data = list(
    K = K,
    N = N,
    T = T,
    N_knots = N_knots,
    N_cores = N_cores,
    y = y,
    rho = rho,
    eta = eta,
    gamma = gamma,
    phi = phi,
    idx_cores = idx_cores,
    d = d,
    d_knots = d_knots,
    d_inter = d_inter,
    w = w,
    P = P, 
    # psi, 
    lag = lag, # not used
    N_p = N_p # not used
  ),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

###### CURRENT ERROR MESSAGES #####
# Running MCMC with 4 parallel chains...
# 
# Chain 1 Exception: mismatch in dimension declared and found in context; processing stage=data initialization; variable name=rho; position=0; dims declared=(11); dims found=(12) (in '/var/folders/v0/n6__xpds389cd0fl9q8yvg3h0000gn/T/RtmpdpIOmx/model-b2945d9a6a0.stan', line 17, column 2 to column 27)
# Warning: Chain 1 finished unexpectedly!
#   
#   Chain 3 Exception: mismatch in dimension declared and found in context; processing stage=data initialization; variable name=rho; position=0; dims declared=(11); dims found=(12) (in '/var/folders/v0/n6__xpds389cd0fl9q8yvg3h0000gn/T/RtmpdpIOmx/model-b2945d9a6a0.stan', line 17, column 2 to column 27)
# Chain 4 Exception: mismatch in dimension declared and found in context; processing stage=data initialization; variable name=rho; position=0; dims declared=(11); dims found=(12) (in '/var/folders/v0/n6__xpds389cd0fl9q8yvg3h0000gn/T/RtmpdpIOmx/model-b2945d9a6a0.stan', line 17, column 2 to column 27)
# Warning: Chain 3 finished unexpectedly!
#   
#   Warning: Chain 4 finished unexpectedly!
#   
#   Chain 2 Exception: mismatch in dimension declared and found in context; processing stage=data initialization; variable name=rho; position=0; dims declared=(11); dims found=(12) (in '/var/folders/v0/n6__xpds389cd0fl9q8yvg3h0000gn/T/RtmpdpIOmx/model-b2945d9a6a0.stan', line 17, column 2 to column 27)
# Warning: Chain 2 finished unexpectedly!
#   
#   Warning: Use read_cmdstan_csv() to read the results of the failed chains.
# Warning messages:
#   1: All chains finished unexpectedly! Use the $output(chain_id) method for more information.
# 
# 2: No chains finished successfully. Unable to retrieve the fit. 

