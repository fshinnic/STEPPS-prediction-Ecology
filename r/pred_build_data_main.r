################################################################################
## This is the starting place for this code repo. It establishes all run parameters
## for one iteration of the STEPPS prediction code and calls the pred_build_data
## code. 
##
## Other Scripts called: Pred_build_data.R
##
## Packages required: None for this script directly.  
##
## Inputs: None. It is just setting run parameters.
##
## Outputs: None specific to this script. All outputs created in Pred_build_data.R


################################IMPORANT ################################################
# CHANGE WHAT LAKES TO RUN ON
lake_group = "wisc_nonvarve_1"  
# lake group options: "wisc_varve", "wisc_nonvarve", "minn_nonvarve"
# single lake optoins: "wisc_nonvarve_1" "wisc_nonvarve_2" "wisc_nonvarve_4"
# minn_nonvarve_1, minn_nonvarve_2, minn_nonvarve_3
################################################################################

# FS- working directory
setwd("/Users/finleyjean/Documents/STEPPS-prediction-Ecology")

#res : relative resolution of grid; should be 1 for 8 km, 3 for 24 km, and 5 for 40 km grids
res  = 3

version="v2.4"
AR     = TRUE #Not sure what this is. In runs some of the 12knots_... files end with AR so maybe autoregressive?
draw   = FALSE #Whether or not to do multiple draws. Code only works if it is TRUE. (pred_build_data 465)
ndraws = 10 #Number of times to run pred_build_data. 
nknots = 120 #Number of knots used in prediction. I do not know if it will work for ones that aren't 120
one_time=FALSE #Whether or not to only predict with one timepoint. (Will break at pred_build_data.r ln 68 if T)
lambda_fixed = TRUE # Is never used
bchron = FALSE #Whether to use bchron or bacon ages 

# FINLEY - change for if dates are varved or not
varve_dates = TRUE

# stat model flags
decomp     = TRUE #Not used
bt         = TRUE #Also never used
mpp        = TRUE
save_plots = TRUE

bacon = TRUE
# draw  = TRUE
add_varves = TRUE
constrain  = FALSE
# how far to extrapolate past last geochron anchor
nbeyond = 5000

age_model  = 'bacon'

suff_dat = '12taxa_mid_comp_ALL_v0.3'


run_pl_Ka_Kgamma = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.4c1', 
                        suff_dat  = suff_dat,
                        kernel    = 'pl', 
                        one_a     = FALSE, 
                        one_b     = TRUE,
                        one_gamma = FALSE, 
                        EPs       = TRUE)

runs = list(run_pl_Ka_Kgamma)

# for (run in runs){
  # for (res in grids){
#     if (draw){
#       for (drN in 1:ndraws){
#         source('r/pred_build_data.r')
#       }
#     }
  # }
# }

# execute all code in the pred_build_data
if (draw){
  for (dr in 1:ndraws){
    source('r/pred_build_data.r')
  }
} else {
  dr = 1
  source('r/pred_build_data.r')
}


