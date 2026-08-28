# instal package to run .stan from Rstudio
setwd("/Users/finleyjean/Documents/STEPPS-prediction-Ecology")

################# NOTES ON DIFFERENCE ON FILE OUTPUT
# what andria did was use rstan - outputs binary file default
# cmdstanr - outputs a .csv
# https://mc-stan.org/cmdstanr/articles/posterior.html
# draws <- posterior::as_draws_rvars(fit$draws())
# x_rvar <- draws$x
# x_array <- posterior::draws_of(draws$x)
####################################################

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


###### LOAD THE ORIGINAL STEPPS DATA PRODUCT
# GOAL - this data product is used to constraint the other data process
# temporal, with both mean and SD at 50 year time intervals
# want to use bayesian model to weight using SD as uncertainty

load("og_STEPPS_products/reference_stepps_grid_cells.RData")
load("og_STEPPS_products/reference_stepps_grid_cells_50yr.RData")

stepps_grid_cells_filtered <- reference_stepps_grid_cells %>% dplyr::filter(
  lake_id == "wisc_nonvarve_2"
) %>% mutate(year = year / 100)

# get STEPPS grid cell in right relative formation
stepps_grid_cells_50yr_filtered <- reference_stepps_grid_cells_50yr %>% dplyr::filter(
  lake_id == "wisc_nonvarve_2"
) %>% 
  dplyr::mutate(year = year / 50) %>% dplyr::select(
  x_center, y_center, lake_id, Taxa, year, reference_prop, reference_sd, taxa_type
)  %>% 
  dplyr::rename("T" =  year) %>% 
  dplyr::mutate(
    T = T-3
  ) %>% 
 dplyr::mutate(
    time_idx = as.integer(T),
    taxa_idx = as.integer(factor(Taxa)),
    cell_idx = as.integer(factor(paste(x_center, y_center)))
  )




####### input parameters #######

# K, N, T, N_cores, N_knots, res, gamma, phi, rho, eta, y, idx_cores, 
# d_knots, d_inter, w, d, lag, P, N_p, meta_pol, meta_pol_all, knot_coords, 
# centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls

# missing the psi parameter, but psi not actually used in the pred_stan_od_mpp_full.stan file

###### compile stan #####

# Compile the model

stan_file <- file.path("cpp/pred_stan_od_mpp_simple.stan")
mod <- cmdstan_model(stan_file) # now works!

############## input files #################

######## MANUALLY CHANGE
#lake_group = "wisc_nonvarve_1"


# lake_groups <- c("minn_nonvarve", "wisc_varve", "wisc_nonvarve", "wisc_nonvarve_1", "wisc_nonvarve_2" ,"wisc_nonvarve_3", "wisc_nonvarve_4" ,  "minn_nonvarve_1", "minn_nonvarve_2", "minn_nonvarve_3")

lake_groups_2 <- c("minn_nonvarve", "wisc_varve", "wisc_nonvarve")
#single_lake_groups <- c("wisc_nonvarve_1", "wisc_nonvarve_2" , "wisc_nonvarve_4" ,
                    #    "minn_nonvarve_1", "minn_nonvarve_2", "minn_nonvarve_3")

for (lake in lake_groups_2) {
  # load in lake group working with:
  if (lake_group == "wisc_varve") {
    load("runs/wisc_varve_1knots_150to2150ybp_PL_PL_PL_PL_test_grid_specs_wisc_varve_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "wisc_nonvarve"){
    load("runs/wisc_nonvarve_4knots_150to2150ybp_PL_PL_test_grid_specs_wisc_nonvarve_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "wisc_nonvarve_1"){
    load("runs/wisc_nonvarve_1_1knots_150to2150ybp_PL_test_grid_specs_wisc_nonvarve_1_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "wisc_nonvarve_2"){
    load("runs/wisc_nonvarve_2_1knots_150to2150ybp_PL_test_grid_specs_wisc_nonvarve_2_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "wisc_nonvarve_3"){
    load("runs/wisc_nonvarve_3_1knots_150to2150ybp_PL_test_grid_specs_wisc_nonvarve_3_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "minn_nonvarve_1"){
    load("runs/minn_nonvarve_1_1knots_150to2150ybp_PL_test_grid_specs_minn_nonvarve_1_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "minn_nonvarve_2"){
      load("runs/minn_nonvarve_2_1knots_150to2150ybp_PL_test_grid_specs_minn_nonvarve_2_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "minn_nonvarve_3"){
    load("runs/minn_nonvarve_3_1knots_150to2150ybp_PL_test_grid_specs_minn_nonvarve_3_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "wisc_nonvarve_4"){
    load("runs/wisc_nonvarve_4_1knots_150to2150ybp_PL_test_grid_specs_wisc_nonvarve_4_v2.4_ar/run1/input.rdata")
  } else if (lake_group == "minn_nonvarve") {
    load("runs/minn_nonvarve_3knots_150to2150ybp_PL_PL_PL_test_grid_specs_minn_nonvarve_v2.4_ar/run1/input.rdata")
  } else {
    stop("Unknown lake_group value")
  }
}


########## RUN STAN ON LOOP FOR ALL SINGLE LAKES################################ 

check_bad <- function(x, name) {
  if (any(is.na(x))) cat(name, "contains NA\n")
  if (any(is.nan(x))) cat(name, "contains NaN\n")
  if (any(is.infinite(x))) cat(name, "contains Inf\n")
}


############## main loop #################

for (lake in lake_groups_2) {


  message("\n==============================")
  message("Running lake group: ", lake)
  message("==============================\n")

  base_path <- paste0(
    "runs/",
    lake,
    "_1knots_150to2150ybp_PL_test_grid_specs_",
    lake,
    "_v2.4_ar/run1/input.rdata"
  )

  # special cases
  if (lake == "minn_nonvarve") {
    base_path <- "runs/minn_nonvarve_3knots_150to2150ybp_PL_test_grid_specs_minn_nonvarve_v2.4_ar/run1/input.rdata"
  }

  if (!file.exists(base_path)) {
    stop("Missing file for lake: ", lake, "\nPath: ", base_path)
  }

  # load files
  load(base_path)

  # preprocessing
  stopifnot(length(rho) == K)
  stopifnot(length(eta) == K)

  rho <- rho[1:(K - 1)]
  eta <- eta[1:(K - 1)]

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

  gamma <- mean(gamma)

  # define weighting matrix
  if (lake %in% c("minn_nonvarve", "wisc_varve")) {
    w <- matrix(1, nrow = 3, ncol = 1)
  } else if (lake == "wisc_nonvarve") {
    w <- matrix(1, nrow = 4, ncol = 1)
  } else if (lake == "minn_nonvarve") {
    w <- matrix(1, nrow = 3, ncol = 1)
  } else if (lake == "wisc_nonvarve_1"){
    w <- matrix(1, nrow = 2, ncol = 1)
  } else {
    w <- matrix(1, nrow = 1, ncol = 1)
  }

  # STEPPS parameters
  stepps_df <- stepps_grid_cells_50yr_filtered

  N_stepps   <- nrow(stepps_df)

  stepps_mean <- stepps_df$reference_prop
  stepps_sd   <- stepps_df$reference_sd

  cell_idx <- stepps_df$cell_idx
  time_idx <- stepps_df$time_idx
  taxa_idx <- stepps_df$taxa_idx

  # run model
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
      lag = lag,
      N_p = N_p,

      # new stepps inputs
      N_stepps = N_stepps,
      stepps_mean = stepps_mean,
      stepps_sd = stepps_sd,
      cell_idx = cell_idx,
      time_idx = time_idx,
      taxa_idx = taxa_idx,
      tau = 3
    ),
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000, #(FINLEY - should change for other ones too)
    output_dir = "/Users/finleyjean/Documents/STEPPS-prediction-Ecology/data/simple_stan_outputs"
  )
  # Warning: 1 of 1 chains had an E-BFMI less than 0.3.
  fit <- fit$summary()


  # save to computer
  output_dir = "/Users/finleyjean/Documents/STEPPS-prediction-Ecology/data/simple_stan_outputs"
  write.csv(
    fit,
    file = file.path(output_dir, paste0("fit_with_og_stpps_", lake, ".csv")),
    row.names = FALSE
  )
}



############ # Run the model for one lake#######################################
  mod_name <- mod$sample(
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

  fit_[lake] <- mod_name$summary()


  setwd("/Users/finleyjean/Documents/STEPPS-prediction-Ecology/data/simple_stan_outputs")
  write.csv(fit_[lake], "fit_", lake, row.names = FALSE)

#
# #### check outputs!!!! #####
#
# # fit_wisc_nonvarve$draws()
# # fit_wisc_nonvarve$draws(format = "df")
# # (fit_wisc_nonvarve$summary())$variable
# # fit_wisc_nonvarve$diagnostic_summary()
#
# fit_wisc_nonvarve_sum <- fit_wisc_nonvarve$summary()
# # head(fit_wisc_nonvarve_sum)
# # Check the column names
# # colnames(fit_wisc_nonvarve_sum)
#
#
# setwd("/Users/finleyjean/Documents/STEPPS-prediction-Ecology/data/simple_stan_outputs")
# write.csv(fit_wisc_nonvarve_sum, "fit_wisc_nonvarve_sum", row.names = FALSE)
#


##### CURRENT STATUS ######
# Rhat (Potential Scale Reduction Factor): Measures convergence by comparing within-chain variance to between-chain variance. A value of Rhat 
# indicates that all chains have converged to the same distribution.
# ess_bulk (Bulk Effective Sample Size): Estimates the effective sample size for the mean and median (center) of the posterior distribution, using rank-normalized draws. It indicates how well the main body of the distribution is sampled.
# ess_tail (Tail Effective Sample Size): Estimates the effective sample size for the variance and 5%/95% quantiles (tails), defined as the minimum of the ESS for both 5% and 95% quantiles. 
# “log density up to a constant”; essentially quantifies how well the model match the data
# Rhats around 2!!!!!!
# lp__
