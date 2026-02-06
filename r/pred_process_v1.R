# Code from r/pred_process_full.r in stepps-prediciton
# appears to be original, un cleaned code for running STEPPS

# all  still exist - 2/6/2026
library(Rcpp)
library(inline)
library(ggplot2)
library(rstan)
library(reshape)
library(fields)

# 
source('r/utils/pred_plot_funs.r')
source('r/utils/pred_helper_funs.r')
source('r/utils/read_stanbin.r')

# suff_dat = '12taxa_457cells_77knots_0to2000ypb_umwE_3by_v0.3'
# suff_fit = '12taxa_457cells_77knots_0to2000ypb_umwE_3by_od_mpp_full'
# suff_dat = '12taxa_459cells_77knots_0to2000ypb_umwW_3by_v0.3'
# suff_fit = '12taxa_459cells_77knots_0to2000ypb_umwW_3by_od_mpp_full'
# suff_dat = '12taxa_459cells_77knots_0to2000ypb_umwW_3by_v0.3'
# suff_fit = '12taxa_459cells_77knots_0to2000ypb_umwW_3by_mpp_full_nug_mu0'
# suff_dat = '12taxa_457cells_77knots_0to2000ypb_umwE_3by_v0.3'
# suff_fit = '12taxa_457cells_77knots_0to2000ypb_umwE_3by_mpp_full_nug_mu0'

# new output scaled res
# suff_dat = 'pred_data_12taxa_387cells_66knots_0to2000ypb_umwE_3by_v0.3'
# suff_fit = '12taxa_387cells_66knots_0to2000ypb_umwE_3by_od_mpp_full_nug_mu0_res'
# 
# suff_fit = '12taxa_459cells_78knots_0to2000ypb_umwW_3by_od_mpp_full_nug_mu0_res'

# suff_dat = '12taxa_699cells_120knots_0to2000ypb_umw_3by_v0.3'
# suff_fit = '12taxa_699cells_120knots_0to2000ypb_umw_3by_od_mpp_full_nug_mu0_res'

# # suff_dat = '12taxa_699cells_120knots_0to2000ypb_PL_umw_3by_v0.3'
# # suff_fit = '12taxa_699cells_120knots_0to2000ypb_PL_umw_3by_tmp'
# 
# suff_dat = '12taxa_699cells_120knots_0to2000ypb_G_umw_3by_v0.3'
# suff_fit = '12taxa_699cells_120knots_0to2000ypb_G_umw_3by_tmp'

# FS chosen files
suff_fit = "120knots_150to2150ybp_PL_test_grid_specs_v2.4_mean_ar"
run = "run1"
load(file.path("runs", suff_fit, run, "input.rdata"))

# where to put the figures
subDir <- paste("figures/", suff_fit, sep='')
suff = ''

# parameters for pred_process_full function
mu0        = TRUE
od         = TRUE
bt         = TRUE
mpp        = TRUE
mut        = FALSE
save_plots = TRUE

# creates path to paste in figures
create_figure_path(subDir)

# FS - need to figure out where this comes from
# Construct file path for the Stan binary output
fname <- sprintf('output/%s.bin', suff_fit)

# Read Stan binary output (custom function)
object <- read_stanbin(fname)

# Extract posterior samples: drop first 4 columns (usually diagnostics) and keep first column at end
samples <- data.frame(object$samples[, 5:ncol(object$samples)], object$samples[, 1])

# Convert to 3D array: iterations x chains x parameters
post <- array(0, c(nrow(samples), 1, ncol(samples)))
post[, 1, ] <- as.matrix(samples)

# Assign parameter names
dimnames(post)[[3]] <- colnames(samples)

# Fallback: if .rdata already exists, load it instead
if (!file.exists(paste0('output/', suff_fit, '.rdata'))) {
  save(post, file = paste0('output/', suff_fit, '.rdata'))
} else {
  load(paste0('output/', suff_fit, '.rdata'))
}


W = K-1

# summary(fit)$summary # take a while

# for full model
N_pars = 3*W + 1

trace_plot_pars(post, N_knots, T, N_pars, taxa=taxa, suff=suff, save_plots=save_plots)
trace_plot_mut(post, N_knots, T, N_pars, mean_type='other', suff=suff, save_plots=save_plots)

## fix this!!!
#trace_plot_knots(fit, N_knots, T, K, N_pars=N_pars, suff=suff, save_plots=save_plots)

# write mean vals to file
# sink(sprintf('%s/%s/summary.txt', wd, path_figs1), type='output')
sink(sprintf('%s/summary.txt', subDir), type='output')
print('The taxa modelled are:')
taxa
cat('\n')
print('Summary of posterior parameter vals:')
get_quants(post, N_pars)
sink()

adj = get_corrections(post, rho, eta, T, K, d_inter, d_knots)
adj_s = adj$adj_s
adj_t = adj$adj_t

####################################################################################################
# chunk: compute and plot proportion chains
####################################################################################################
process_out = build_props_full(post, rho, eta, T, K, d, d_inter, d_knots, od, mpp, mu0)
# save(process_out, file='r/pred/dump/process_out.rdata')
# rm(post)


# M = diag(N*T)
# M = M - P
# 
# foo = M %*% Halpha[,1,1]
# 
# # Halpha

mu_g   = process_out$mu_g
r_pred = process_out$r
g      = process_out$g
Halpha_t = process_out$Halpha_t
Halpha_s = process_out$Halpha_s
# rm(process_out)

niter = dim(post[,1,])[1]
mean_Halpha_t = array(NA, dim=c(W, T-1, niter))

for (k in 1:W){
  for (t in 1:(T-1)){
    mean_Halpha_t[k, t, ] = colSums(Halpha_t[,(k-1)*(T-1)+t,])/N
  }
}

mu_t = get_mut(post, N_pars, W)

pdf(file=paste0(subDir, '/trace_Halpha_t.pdf'), width=8, height=12)
for (k in 1:W){
  par(mfrow=c(5,2))
  par(oma=c(0,0,2,0))
  for (t in 1:(T-1)){
    plot(mean_Halpha_t[k,t,], type="l", ylab=paste0('mean_Halpha_t[', k, ',', t, ']'))
    abline(h=mean(mean_Halpha_t[k,t,]), col="blue")
    abline(h=quantile(mean_Halpha_t[k,t,],probs=0.025), col='blue', lty=2)
    abline(h=quantile(mean_Halpha_t[k,t,],probs=0.975), col='blue', lty=2)
    #   lines(mu_t[t+1,1,], col="blue")
    #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
    
  }
  #     title(main=taxa[k], outer=TRUE)
}
dev.off()


niter = dim(post[,1,])[1]
mean_Halpha_s = array(NA, dim=c(W, niter))

for (k in 1:W){
  mean_Halpha_s[k, ] = colSums(Halpha_s[,k,])/N
}

mu = get_mu(post, W)
pdf(file=paste0(subDir, '/trace_mu.pdf'), width=8, height=12)
par(mfrow=c(5,2))
par(oma=c(0,0,2,0))
for (k in 1:W){
  plot(mu[k,], type="l", ylab=paste0('mu[', k, ']'))
  abline(h=mean(mu[k,]), col="blue")
  abline(h=quantile(mu[k,],probs=0.025), col='blue', lty=2)
  abline(h=quantile(mu[k,],probs=0.975), col='blue', lty=2)
  #   lines(mu_t[t+1,1,], col="blue")
  #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
  #     title(main=taxa[k], outer=TRUE)
}
dev.off()

pdf(file=paste0(subDir, '/trace_Halpha_s.pdf'), width=8, height=12)
par(mfrow=c(5,2))
par(oma=c(0,0,2,0))
for (k in 1:W){
  plot(mean_Halpha_s[k,], type="l", ylab=paste0('mean_Halpha_s[', k, ']'))
  abline(h=mean(mean_Halpha_s[k,]), col="blue")
  abline(h=quantile(mean_Halpha_s[k,],probs=0.025), col='blue', lty=2)
  abline(h=quantile(mean_Halpha_s[k,],probs=0.975), col='blue', lty=2)
  #   lines(mu_t[t+1,1,], col="blue")
  #     plot(sum_Halpha_t[k,t,], type="l", ylab=paste0('sum_Halpha[', k, ',', t, ']'))
  #     title(main=taxa[k], outer=TRUE)
}
dev.off()


# ######################################
# # check if need further OD
# mu_t       = get_mut(post, N_pars=W)
# sum_Halpha = process_out$sumHalpha
# 
# pdf(file=paste0(subDir, '/compare_mu_Halpha.pdf'), width=8, height=8)
# for (k in 1:W){
#   par(mfrow=c(3,2))
#   par(oma=c(0,0,2,0))
#   for (t in 1:3){
#     plot(mu_t[t,k,], type="l", ylab=paste0('mu_t[', t, ',', k, ']'))
#   #   lines(mu_t[t+1,1,], col="blue")
#     plot(sum_Halpha[t,k,], type="l", ylab=paste0('sum_Halpha[', t, ',', k, ']'))
# 
#   }
#   title(main=taxa[k], outer=TRUE)
# }
# dev.off()
# 
# mu_t_int = colSums(mu_t[,1,])
# 
# mu = get_mu(post, N_pars=W)
# 
# mu_t_int - mu[,1]
# 
# pdf(file=paste0(subDir, '/compare_mus.pdf'), width=8, height=6)
# for (k in 1:W){
#   mu_t_int = colSums(mu_t[,k,])
#   par(mfrow=c(2,1))
#   plot(mu_t_int, type='l', ylab=paste0('int(mu_t[t, ', k, '])'))
#   plot(mu[,k], type='l',  ylab=paste0('mu[', k, ']'), col='black')
# }
# dev.off()
#######################################

trace_plot_process(mu_g, suff='mu_g', save_plots=save_plots)
trace_plot_process(r_pred, suff='r', save_plots=save_plots)
trace_plot_process(g, suff='g', save_plots=save_plots)


r_mean = matrix(NA, nrow=N*T, ncol=K)
g_mean = matrix(NA, nrow=N*T, ncol=W)
niter = dim(r_pred)[3]

for (i in 1:(N*T)){
  r_mean[i,] = rowSums(r_pred[i,,])/niter
  g_mean[i,] = rowSums(g[i,,])/niter
}

rm(g)
rm(r_pred)


limits = get_limits(centers_pls)
####################################################################################################
# chunk: plot predicted distributions
####################################################################################################
suff1=paste(suff_fit, '_props', sep='')

# plot_pred_maps(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1,  save_plots=save_plots)

suff1.1=paste(suff_fit, '_props_select', sep='')
plot_pred_maps_select(r_mean, centers_veg, taxa=taxa, ages, N, K, T, thresh=0.5, limits, type='prop', suff=suff1.1,  save_plots=save_plots)

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# p_binned <- plot_pred_maps_binned(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff=suff1, save_plots=save_plots)
p_binned <- pred_maps_binned_select(r_mean, centers_veg, breaks, taxa, ages, N, K, T, limits, suff1, save_plots, fpath=subDir)
####################################################################################################
# chunk: predicted process maps
####################################################################################################
suff2=paste(suff_fit, '_process', sep='')

plot_pred_maps(g_mean, centers_veg, taxa=taxa, ages, N, K-1, T, thresh=NA, limits, type='process', suff=suff2,  save_plots=save_plots)

####################################################################################################
# chunk: plot observed proportions
####################################################################################################
suff3=paste(suff_fit, '_data', sep='')

plot_data_maps(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, thresh=0.5, limits, suff=suff3, save_plots=save_plots)

suff4=paste(suff_fit, '_data_binned', sep='')
plot_data_maps_binned(y_veg, centers=centers_pls, taxa=taxa, ages, N_pls, K, T, breaks, limits, suff=suff4, save_plots=save_plots)


# N_cores and centers_polU broken for split domain data....
idx.keep  = c(1,2,length(ages)/2,T)
# plot_core_locations(y, centers_polU, centers_pls, ages, limits)
plot_core_locations_select(y, centers_pol, centers_pls, ages, idx.keep, limits, fpat=subDir)


# plot5<-grid.arrange(arrangeGrob(c, p_binned, nrow=1), 
#                      widths=c(10))
# plot5<-grid.arrange(c, p_binned, heights=c(10,10), nrow=1)


# 
# # objects are c, d, p_binned
# # Get the widths
# g <- ggplotGrob(p_binned)
# keep <- !grepl("strip-top", g$layout$name)
# g$grobs <- g$grobs[keep]
# g$layout <- g$layout[keep, ]
# 
# # gA <- ggplot_gtable(ggplot_build(p_binned))
# # gB <- ggplot_gtable(ggplot_build(c))
# 
# c_grob <- ggplotGrob(c)
# c_grob$heights <- g$heights
# 
# # gtable_show_layout(g)
# gf <- gtable_add_cols(g, unit(5, "cm"), pos=0)
# gf <- gtable_add_cols(g, unit(0.25, "lines"), pos=0)
# gf <- gtable_add_grob(gf, c_grob,
#                      t = 3, l=1, b=nrow(g), r=1)
# g <- gtable_add_rows(g, unit(8, "cm"), pos=0)
# # g <- gtable_add_rows(g, g$heights[1], pos=0)
# g <- gtable_add_grob(g, ggplotGrob(d),
#                      t = 1, l=2, b=1, r=ncol(g)-2)
# grid.newpage()
# pdf('figures/pred/test_gtable.pdf', width=12, height=16)
# grid.draw(g)
# dev.off()
# # 
# 
# 
# ng = nullGrob()
# g <- gtable(widths = unit(c(5, 7), "in"),  # need lcm(3,2)=6 for the matrix rows
#             heights = unit(c(5, 7), "in"))
# gtable_show_layout(g)
# 
# # # Add a grob:
# # rect <- rectGrob(gp = gpar(fill = "black"))
# # a <- gtable_add_grob(a, rect, 1, 1)
# # a
# # plot(a)
# 
# g<- gtable_add_grob(g, c_grob, 2, 1)
# g<- gtable_add_grob(g, ggplotGrob(d), 1, 2)
# g<- gtable_add_grob(g, ggplotGrob(p_binned), 2, 2)
# 
# grid.newpage()
# grid.draw(g)


# x1 <- freeze_panels(d, draw=TRUE, width=unit(1,"in"),height=unit(1,"in"))
# x2 <- freeze_panels(c, draw=TRUE, width=unit(1,"in"),height=unit(1,"in"))
# x3 <- freeze_panels(p_binned, draw=TRUE, width=unit(1,"in"),height=unit(1,"in"))
# 
# 
# 
# 
# 
# grid.arrange(arrangeGrob(ng, x1, nrow=1, widths=c(0.5,3.5)),
#              arrangeGrob(x2, x3, nrow=1, widths=c(0.5,3.5)),
#              nrow=2)
# 
# grid.arrange(ng, x1, x2, x3, nrow=2, widths=c(0.5,3.5))
# 
# 
# library(grid)
# ht = 3
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(ht, length(ages*ht)), "null"))))
# grid.rect()
# # grid.text("title of this panel", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
# print(d, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(c, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(p_binned, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
# 
# w = 5
# grid.newpage()
# pushViewport(viewport(width=10, height=10, default.units="in")); 
# #gt = grobTree(rectGrob(), circleGrob()) 
# grid.arrange(c, p_binned, nrow=1, newpage=FALSE,
#              widths = unit(c(w, w*7), "null"))
# upViewport()
# 
# 
# top <- arrangeGrob(r, r, r, r, nrow=2)
# 
# grid.arrange(top, r, ncol=1, main = "Main Title", heights=c(2, 1))
# 
# grid.newpage()
# align_plots(c,p_binned)
# # ####################################################################################################
# # # chunk: plot predicted and observed proportions in same frame
# # ####################################################################################################
# # suff4=paste(suff, '_compare', sep='')
# # 
# # plot_both_maps(r_mean, y_veg, centers=centers_pls, taxa=taxa, ages, N, K, T, thresh=1, suff=suff3, save_plots=save_plots)
