# Code from r/pred_process_full.r in stepps-prediciton
# appears to be original, un cleaned code for running STEPPS

# Original from pred_process_full.r
#process_out = build_props_full(cal_post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0)
# edited to match file times

mu0        = TRUE # required in function, made true from old code
process_out = build_props_full(cal_post, rho, eta, T, K, d, d_inter, d_knots, mpp, mu0)
