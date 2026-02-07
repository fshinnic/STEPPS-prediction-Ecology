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
suff_fit = "120knots_150to2150ybp_PL_test_grid_specs_v2.4_mean_ar"
run = "run1"
load(file.path("runs", suff_fit, run, "input.rdata"))

####### input parameters #######
#K, N, T, N_cores, N_knots, res, gamma, phi, rho, eta, y, idx_cores, d_knots,
#d_inter, w, lag, P, N_p, meta_pol, meta_pol_all, knot_coords, centers_pls, 
# centers_veg, centers_pol, taxa, ages, y_veg, N_pls

