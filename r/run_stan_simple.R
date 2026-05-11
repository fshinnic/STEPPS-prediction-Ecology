# instal package to run .stan from Rstudio
setwd("/Users/finleyjean/Documents/STEPPS-prediction-Ecology")
library(cmdstanr)
#cmdstanr::install_cmdstan()

####### Old GOAL ####### 
# Run STAN using one core (LilyLake) and one grid cell
# Hopefully reduce the spatial importance and see if i can just use the code to get posteriors for vegetation over time

####### Current GOAL ####### 

# Run STAN using 3 cors from one grid cell to have more pollen information (check if improves rhat for convergence)

####### END GOALs ####### 

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

############## input files #################

# FS chosen files (created in pred_build_data_main.r calling of pred_build_data.r)
# fname : "runs/1knots_150to2150ybp_PL_test_grid_specs_v2.4_ar/run1/input"
# dirName = "runs/1knots_150to2150ybp_PL_test_grid_specs_v2.4_ar"
# subDir = "run1"
# fname <- file.path(dirName, subDir, 'input')  # same as used when saving
# rdata_file <- paste0(fname, '.rdata')
# 
# # Load the file
# load(rdata_file)

# MANUALLY CHANGE
lake_group = "wisc_varve"
# options: "minn_nonvarve", "wisc_varve", "wisc_nonvarve"

# load in lake group working with:
if (lake_group == "wisc_varve") {
  load("runs/wisc_varve_1knots_150to2150ybp_PL_test_grid_specs_wisc_varve_v2.4_ar/run1/input.rdata")
  
} else if (lake_group == "wisc_nonvarve") {
  load("runs/wisc_nonvarve_4knots_150to2150ybp_PL_test_grid_specs_wisc_nonvarve_v2.4_ar/run1/input.rdata")
  
} else if (lake_group == "minn_nonvarve") {
  load("runs/minn_nonvarve_3knots_150to2150ybp_PL_test_grid_specs_minn_nonvarve_v2.4_ar/run1/input.rdata")
  
} else {
  stop("Unknown lake_group value")
}
####### input parameters #######

# K, N, T, N_cores, N_knots, res, gamma, phi, rho, eta, y, idx_cores, 
# d_knots, d_inter, w, d, lag, P, N_p, meta_pol, meta_pol_all, knot_coords, 
# centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls

# missing the psi parameter, but psi not actually used in the pred_stan_od_mpp_full.stan file

###### compile stan #####

# Compile the model
stan_file <- file.path("cpp/pred_stan_od_mpp_simple.stan")
mod <- cmdstan_model(stan_file) # now works!


######### fix 1
# The stan code appears to have a reference taxa
# takes rho = k-1, where right now the input params. have rho = k (12)
# drop the "other" taxa to serve as a reference

stopifnot(length(rho) == K)
stopifnot(length(eta) == K)

rho <- rho[1:(K-1)]
eta <- eta[1:(K-1)]


########## fix 2
# check if the data contains a NA or missing value or an infemum
check_bad <- function(x, name) {
  if (any(is.na(x)))  cat(name, "contains NA\n")
  if (any(is.nan(x))) cat(name, "contains NaN\n")
  if (any(is.infinite(x))) cat(name, "contains Inf\n")
}

check_bad(rho, "rho")
check_bad(eta, "eta")
check_bad(phi, "phi")
check_bad(y, "y")
check_bad(d, "d")
check_bad(d_knots, "d_knots")
check_bad(d_inter, "d_inter")
check_bad(w, "w")
check_bad(P, "P")
check_bad(gamma, "gamma")

gamma = mean(gamma)

# issue seems to be that the original stsan model expected gamma = scalar
# changed gamma to be vector, calling gamme[K] since taxa-specific

########## fix 3
# variable name = w
# dims declared = (199, 6378)

# choice, remove 12 w = (12, 199, 6378)

# FINLEY - Come back to this
# reduces information by taxa (declared (1,17), found (12,1, 127)[...])
# w <- matrix(1, nrow = 1, ncol = 1)###### getting non-posistve definite matrix , having issues with it sampeling across space
w <- matrix(1, nrow = 3, ncol = 1) 

# print(rho)
# summary(rho)        # check for very small or huge values
# summary(d_knots)    # check for zeros or Inf/NaN
w

########## RUN STAN 
# Run the model
fit_w_one_v4 <- mod$sample(
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
  parallel_chains = 1,
  iter_warmup = 1000,
  iter_sampling = 10,
  output_dir = "/Users/finleyjean/Documents/STEPPS-prediction-Ecology/data/simple_stan_outputs"
)

#### check outputs!!!! #####

fit_w_one_v4$draws()
fit_w_one_v4$draws(format = "df")
(fit_w_one_v4$summary())$variable
fit_w_one_v4$diagnostic_summary()

fit_w_one_v4 <- fit_w_one_v4$summary()  
head(fit_w_one_v4)
# Check the column names
colnames(fit_w_one_v4)

View(fit_w_one_v4)

setwd("/Users/finleyjean/Documents/STEPPS-prediction-Ecology/data/simple_stan_outputs")
write.csv(fit_w_one_v4, "stan_simple_varve_v4", row.names = FALSE)



##### CURRENT STATUS ######
# Rhat (Potential Scale Reduction Factor): Measures convergence by comparing within-chain variance to between-chain variance. A value of Rhat 
# indicates that all chains have converged to the same distribution.
# ess_bulk (Bulk Effective Sample Size): Estimates the effective sample size for the mean and median (center) of the posterior distribution, using rank-normalized draws. It indicates how well the main body of the distribution is sampled.
# ess_tail (Tail Effective Sample Size): Estimates the effective sample size for the variance and 5%/95% quantiles (tails), defined as the minimum of the ESS for both 5% and 95% quantiles. 
# “log density up to a constant”; essentially quantifies how well the model match the data
# Rhats around 2!!!!!!
# lp__
