library(fields)
library(rstan)
library(sp)
#library(rgdal)
library(ggplot2)
library(mvtnorm)
#library(maptools)
library(maps)
library(plyr)
library(sf)


source(file.path(getwd(), "r/utils/pred_helper_funs.r"))

getwd()

######################################################################################################################################
# user defs
######################################################################################################################################

# us albers shape file
#us.shp <- readOGR('data/map/us_alb.shp', 'us_alb')
us.shp <- sf::st_read('data/map/us/us_alb.shp', 'us_alb')

# reconstruction limits and bin-width
int  = 100
tmin = 150
if (one_time) {
  tmin = tmin + (slice-1)*int
  tmax = tmin + int
} else {
  tmax = tmin + 20*int
}

# rescale
rescale = 1e6

# suffix to append to figure names
#HO - I have no idea where grid_specs comes from so I will add my own to get it to run.
grid_specs = "test_grid_specs"
suff    = paste(grid_specs, '_', version, sep='') 

# specify states; relict from testing
states_pol = c('minnesota', 'wisconsin', 'michigan:north')
states_pls = c('minnesota', 'wisconsin', 'michigan:north')

# specify the taxa to use
# must be from the list: taxa_sub = c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beeh', 'elm', 'spruce', 'ash', 'hemlock')
# always have: 'other.hardwood' and 'other.conifer'
taxa_all = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock'))
taxa_sub = toupper(c('oak', 'pine', 'maple', 'birch', 'tamarack', 'beech', 'elm', 'spruce', 'ash', 'hemlock'))

# K is number of taxa
K = as.integer(length(taxa_sub) + 1)
W = K-1

##########################################################################################################################
## paths and filenames to write to meta file
##########################################################################################################################
suff_veg = paste0('12taxa_6341cells_', nknots, 'knots')

path_grid = 'data/grid/umw_3by_v2.rdata'

path_pls      = 'data/pls_umw_v0.6.csv'

path_pollen   = 'data/sediment_ages_v1.0_varves.csv'

# Only non-varve lakes of interest
# path_pollen = 'data/fs_data/wisc_nonvarve_pollen_ts.csv'

if (bchron){
  path_age_samples    = 'data/bchron_ages'
} else {
  path_age_samples    = 'data/bacon_ages'
}

path_cal      = paste0('data/calibration/', runs[[1]]$suff_fit,'.csv')
path_veg_data = paste0('data/veg_data_', suff_veg, '_v0.4.rdata')


# path_veg_pars = paste0('data/', suff_veg, '_nb_v0.5/veg_pars_', nknots, 'knots.rdata') #HO -This is not a file or directory that exists. Replacing with what I think is correct
path_veg_pars = 'data/veg_pars_120knots.rdata'
##########################################################################################################################
## read in tables and data
##########################################################################################################################

# conversion tables
tree_type = read.table('data/assign_HW_CON.csv', sep=',', row.names=1, header=TRUE)
convert   = read.table('data/dict-comp2stepps.csv', sep=',', row.names=1, header=TRUE)

# pls veg data
pls.raw = data.frame(read.table(file=path_pls, sep=",", row.names=1, header=TRUE))

# read in grid
#load(file=path_grid)

# pollen data
pollen_ts = read.table(path_pollen, header=TRUE, sep=',', stringsAsFactors=FALSE)

# FS - each lake has unique id and multiple bacon_draw#
pol_ids = data.frame(id=unique(pollen_ts$id), stat_id=seq(1, length(unique(pollen_ts$id))))

# if draw=TRUE then replace mean age_bacon with draw age_bacon
if (draw) {
  all_files   = list.files(path_age_samples)
  all_drawRDS = all_files[grep('draw', all_files)]
  
  # FS - randomply choses one draw index
  drawRDS     = all_drawRDS[sample(seq(1, length(all_drawRDS)), 1)]   # random sample from available posterior age draws
  
  age_sample   = readRDS(file.path(path_age_samples, drawRDS))
  # replace age_bacon with the draw
  if (bchron){
    pollen_ts$age_bchron = age_sample
  } else {
    pollen_ts$age_bacon  = age_sample  
  }
}

if (bchron){
  pollen_ts = pollen_ts[!is.na(pollen_ts$age_bchron),]
} else {
  pollen_ts = pollen_ts[!is.na(pollen_ts$age_bacon),]
}


# max_ages
if (constrain){
  if (bchron){
    max_ages = read.table(file='data/pol_ages_bchron_6.csv', sep=',', header=TRUE)
  } else {
    max_ages = read.table(file='data/pol_ages_v6.csv', sep=',', header=TRUE)
  }
  drop_samples = constrain_pollen(pollen_ts, max_ages, nbeyond=nbeyond)
  
  if (add_varves){
    vids = c(2309, 14839, 3131)
    drop_samples[which(pollen_ts$id %in% vids)] = FALSE
  }
  pollen_ts = pollen_ts[!drop_samples,]
}


cal_fit = rstan::read_stan_csv(path_cal)

# read in veg data and output
# veg data specifies which knots to use
load(path_veg_data)
veg_post = readRDS(file=path_veg_pars)
# load(file=path_veg_pars)
# veg_post = post

##########################################################################################################################
## read in and organize pls data
##########################################################################################################################

# changes taxa names to lower case
colnames(pls.raw) = tolower(colnames(pls.raw))

# pull the subset of proportions
taxa.start.col = min(match(tolower(rownames(convert)), colnames(pls.raw)), na.rm=TRUE)

# # might need to fix this later, doesn't work with updated data but don't need it now
# if (any(!(tolower(sort(taxa)) == sort(colnames(pls_dat))))) {
#   pls_dat  = pls.raw[,taxa.start.col:ncol(pls.raw)]
#   colnames(pls_dat) = as.vector(convert[match(colnames(pls_dat), tolower(rownames(convert))),1])
#   pls_dat_collapse  = sapply(unique(colnames(pls_dat)), 
#                              function(x) rowSums( pls_dat[ , grep(x, names(pls_dat)), drop=FALSE]) )
#   counts = data.frame(pls_dat_collapse[,sort(colnames(pls_dat_collapse))])
# }

# FS - raw pollen couts for each taxa
counts = pls.raw[,taxa.start.col:ncol(pls.raw)]

# FS - geographical location of each count
meta   = pls.raw[,1:(taxa.start.col-1)]
# kilometers
# pls$X = pls$X/1000
# pls$Y = pls$Y/1000

# FS - use the new function (inputed to pred_hilders_funs.r)
# create collumns state (minnesota, south dakatoa, north dakota, iowa, na, wiconsin, michigan:north, illinois, michigan:south, indiana, ohio)
# and state2 (minnesota, wisconsin, michigan:north, and michigan:south) the corrected map look up
meta        = split_mi(meta)

# only keep counts/meta that are in the corrected "minnesota" "wisconsin" "michigan:north" (states_pls)
counts      = counts[which(meta$state2 %in% states_pls),]
meta        = meta[which(meta$state2 %in% states_pls),]


centers_pls = data.frame(x=meta$x, y=meta$y)/rescale # megameters!

# FS - Corrected for SP object; checks that seperated out the correct uppermidwest cells
# plot(centers_pls[,1]*rescale, centers_pls[,2]*rescale, asp=1, axes=F,  col='antiquewhite4', xlab='',ylab='', pch=19, cex=0.2) # old code
# plot(us.shp, add=T) # old code
# plots center of pls grid cells (8 km x 8 km)
plot(centers_pls[,1] * rescale,
     centers_pls[,2] * rescale,
     asp = 1,
     axes = FALSE,
     col = "antiquewhite4",
     pch = 19,
     cex = 0.2)
# adds borders
plot(st_geometry(us.shp), add = TRUE, border = "black", lwd = 0.5)

# FS - aggreg raw pollen count data into 12 taxa (including OTHER group)
y_veg = convert_counts(counts, tree_type, taxa_sub)

# FS - removes collumn names, creating counts as ordered matrix
taxa = colnames(y_veg)
y_veg = as.matrix(round(unname(y_veg)))
rownames(y_veg) = NULL
y_veg = unname(y_veg)

# y = y_build(counts, taxa_sub) # fix this if we want to use a subset of taxa

K = as.integer(ncol(y_veg))
W = K-1
N_pls = nrow(y_veg)

# make sure columns are in order! 
# y_veg = y_veg[,taxa]

##########################################################################################################################
## chunk: read in coarse grid and pollen data
##########################################################################################################################
# ---- FS made domain wisconsin? 
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthdata)
# 
# # Get US states
# us_states <- ne_states(country = "United States of America", returnclass = "sf")
# 
# # Subset Wisconsin
# wi <- us_states[us_states$name == "Wisconsin", ]
# 
# # Transform to Albers Equal Area (EPSG 5070)
# wi_albers <- st_transform(wi, 5070)
# 
# # Extract coordinates as a data frame
# domain <- st_coordinates(wi_albers)[, 1:2]
# domain <- data.frame(lat = domain[,1], long = domain[,2])
# 
# head(domain)
### ---- end of FS work

# FS - created domain of uppermidwest form PLS grid centers
states_pls <- c("wisconsin", "michigan:north","minnesota")
domain <- meta[meta$state2 %in% states_pls, c("x", "y")]

# FIXME: ADD STATE TO GRID
# coarse_domain  = coarse_domain[coarse_domain$state %in% states_pls,]
coarse_centers = domain[,1:2]

# check domain working with
plot(coarse_centers[,1] * rescale,
     coarse_centers[,2] * rescale,
     col = "blue",
     pch = 19,
     cex = 0.3,
     asp = 1,
     axes = FALSE,
     xlab = "",
     ylab = "")
plot(st_geometry(us.shp), add = TRUE, border = "black")



# assign grid to centers_veg
centers_veg = coarse_centers
N = nrow(centers_veg)

# subdomain boundaries
xlo = min(centers_veg$x)
xhi = max(centers_veg$x)
ylo = min(centers_veg$y)
yhi = max(centers_veg$y)

##########################################################################################################################
## chunk: reorganize pollen data
##########################################################################################################################

# set tamarack to 0 at tamarack creek; see Dawson et al. QSR 2016
pollen_ts[pollen_ts$id == 2624, 'TAMARACK'] = rep(0, sum(pollen_ts$id == 2624))

saveRDS(pollen_ts, file='data/pollen_ts.RDS')

# FS - only keep sites in "minnesota"      "wisconsin"      "michigan:north"
pollen_ts1 = pollen_ts[which(pollen_ts$state %in% states_pol),]
colnames(pollen_ts1)

# FS - had to edit this function to use sf instead (rerun pollen_to_albers.r file)
# reproject pollen coords from lat long to Albers
pollen_ts2 = pollen_to_albers(pollen_ts1)
# FS - location of pollen
pollen_locs = cbind(pollen_ts2$x, pollen_ts2$y)

# FS - checks that the pollen sitees are within the vegetation domain
plot(centers_veg$x, centers_veg$y, col='green', pch=19, main='Vegetation vs Pollen')
points(pollen_locs[,1], pollen_locs[,2], col='red', pch=19)


# FS - FIX ME - currently no matches ie no cores in the domain
#pollen_int  = cores_near_domain(pollen_locs, centers_veg, cell_width = res*8000/rescale)
cell_width_value <- 8000   #  FS - changed!!!!
pollen_int <- cores_near_domain(pollen_locs, centers_veg, cell_width = cell_width_value)

# FS - Keep only pollen points that exactly match the coordinates in pollen_int
# FS - idx_pollen_int is a logical vector used to filter pollen_ts2 to these points
idx_pollen_int = apply(pollen_locs, 1, 
                       function(x) if (any(rdist(x, pollen_int) < 1e-8)) {return(TRUE)} else {return(FALSE)})
pollen_ts3 = pollen_ts2[idx_pollen_int, ]

# FS - plot domain and core locations in Wisconsin, upper michican, and minnesota
par(mfrow=c(1,1))
plot(centers_veg$x*rescale, centers_veg$y*rescale)
points(pollen_ts3$x*rescale, pollen_ts3$y*rescale, col='blue', pch=19)
plot(us.shp, add=T, lwd=2)

##########################################################################################################################
## chunk: prepare pollen data; aggregate over time intervals
##########################################################################################################################
# FS - Moved x and y columns to front so that they aren't lost in build_pollen_counts function
pollen_ts3 <- pollen_ts3[, c("id", "x", "y", setdiff(colnames(pollen_ts3), c("id","x","y")))]

# sum counts over time interval length intervals
pollen_agg = build_pollen_counts(tmin=tmin, tmax=tmax, int=int, pollen_ts=pollen_ts3, taxa_all, taxa_sub, age_model=age_model)
#pollen_agg = build_pollen_counts_fast_core(tmin=tmin, tmax=tmax, int=int, pollen_ts=pollen_ts)

# saveRDS(pollen_ts3, file=paste0(subDir, '/pollen_meta.RDS'))

meta_pol_all = pollen_agg[[3]]
meta_pol   = pollen_agg[[2]]
counts     = pollen_agg[[1]]

meta_pol$stat_id = pol_ids$stat_id[match(meta_pol$id, pol_ids$id)]
meta_pol_all$stat_id = pol_ids$stat_id[match(meta_pol_all$id, pol_ids$id)]

pollen_ts$stat_id = pol_ids$stat[match(pollen_ts$id, pol_ids$id)]

ages    = unique(sort(meta_pol$age)) # FS - only 21, each = 100 years

# number of time slaces (20 of 100 years = 2000 years)
T       = length(ages) 

if (one_time) {
  lag = 0
} else {
  lag     = unname(as.matrix(dist(matrix(ages), upper=TRUE)))
}

N_cores = length(unique(meta_pol$id)) # FS - 199

y = convert_counts(counts, tree_type, taxa_sub)

# make sure columns match!
if (sum(colnames(y) %in% taxa) != K){
  print('The number of taxa wanted does not match the number of taxa in the data frame! Name mismatch likely.')
}

# y = y[,taxa]
y = unname(y)

centers_pol = data.frame(x=numeric(N_cores), y=numeric(N_cores))

for (i in 1:N_cores){
  id = unique(meta_pol$id)[i]
  idx = min(which(meta_pol$id == id))
  print(idx)
  centers_pol[i,] = c(meta_pol$x[idx], meta_pol$y[idx])  
}


# some are duplicates, but we still need them as separate rows!
# centers_pol <- meta_pol[!duplicated(cbind(meta_pol$x, meta_pol$y)), c('x', 'y')]

# indices for which cells the cores fall in
idx_cores <- build_idx_cores(centers_pol, centers_veg, N_cores)

plot(centers_veg$x*rescale, centers_veg$y*rescale, col='lightgrey')
points(centers_veg[idx_cores,'x']*rescale, centers_veg[idx_cores,'y']*rescale, col='red', pch=19)
points(centers_pol$x*rescale, centers_pol$y*rescale, col='blue', pch=4, cex=1.4)
#plot(us.shp, add=TRUE)

## SAVE ENVIORNMENT THERE TO NOT RUN POLLEN AGG AGAIN ##
# save.image("my_environment.RData")
# load("my_environment.RData")

# check domain splitting
# FS - pollen check doesn't exist
# idx_cores_all <- build_idx_cores(cbind(pollen_check$x, pollen_check$y), centers_veg, N_cores=nrow(pollen_check))

c##########################################################################################################################
## chunk 3: build distance matrices
##########################################################################################################################

# matrix between all possible begetation squares
d = rdist(centers_veg, centers_veg)
diag(d) <- 0

d_knots = rdist(knot_coords, knot_coords)
diag(d_knots) <- 0

d_inter = rdist(centers_veg, knot_coords)
d_inter[which(d_inter<1e-8)]=0

d_pol = rdist(centers_pol, centers_veg)
d_pol[which(d_pol<1e-8)]=0

N_knots = nrow(knot_coords)

##########################################################################################################################
## pull in calibration parameters
##########################################################################################################################

KW     = FALSE
KGAMMA = FALSE

# kernel    = run$kernel
# FS -  to check that there is only one kernel object
kernel <- sapply(runs, function(x) x$kernel)

# phi = differential production
# gamma = proportion of that cores pollen the given taxa producted in the grid cell
# Log_a, mu_gamma, sigma_gamma,b, a, not
cal_post      = rstan::extract(cal_fit, permuted=FALSE, inc_warmup=FALSE)
col_names = colnames(cal_post[,1,])
par_names  = unlist(lapply(col_names, function(x) strsplit(x, "\\[")[[1]][1]))

if (draw) {
  draw_cal = sample(seq(1, dim(cal_post)[1]), 1)
  cal_post     = cal_post[draw_cal,1,]
} else {
  cal_post = colMeans(cal_post[,1,])
}

phi = unname(cal_post[which(par_names == 'phi')][1:K])

one_gamma = sapply(runs, function(x) x$one_gamma) 

if (one_gamma){
  gamma = unname(cal_post[which(par_names == 'gamma')])
} else {
  KGAMMA = TRUE
  gamma  = unname(cal_post[which(par_names == 'gamma')][1:K])
}

if (kernel=='pl'){
  
  # FS
  # one_a = runs$one_a
  one_a = sapply(runs, function(x) x$one_a) 
  
  if (one_a){
    a = unname(cal_post[which(par_names == 'a')])
  } else {
    KW = TRUE
    a = unname(cal_post[which(par_names == 'a')][1:K])
  }
  
  # FS
  # one_b = runs$one_b
  one_b = sapply(runs, function(x) x$one_b) 
  
  if (one_b){
    b = unname(cal_post[which(par_names == 'b')])
  } else {
    KW = TRUE
    b = unname(cal_post[which(par_names == 'b')][1:K])
  }
}

# FS - START HERE 1/15/2025 - Check more what this function does

# FS - changed since runs appears to have multiple layers
# create the weighting matrix for pollen dispersion for each of the 12 taxa groups
w <- build_weight_matrix(cal_post, d_pol, idx_cores, N, N_cores, runs[[1]])

#####################################################################################
# calculate potential d
# used to determine C normalizing constant in the non-local contribution term
#####################################################################################

coord_pot = seq(-700000, 700000, by=8000)
coord_pot = expand.grid(coord_pot, coord_pot)

d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)
d_pot = unname(as.matrix(count(data.frame(d_pot))))

N_pot     = nrow(d_pot)

#####################################################################################
# recompute gamma; needed to account for change in resolution from base res
#####################################################################################

# FS - they made the resolution coarses, but d_hood doesn't exist
# w_coarse  = build_sumw_pot(cal_post, K, length(d_hood), cbind(t(d_hood), rep(1, length(d_hood))), run)
# gamma_new = recompute_gamma(w_coarse, sum_w_pot, gamma)

#####################################################################################
# veg run pars
#####################################################################################
par_names = sapply(strsplit(colnames(veg_post), '\\.'), function(x) x[[1]])
eta = veg_post[,which(par_names == 'eta')]
rho = veg_post[,which(par_names == 'rho')]

if (draw){
  iter = sample(seq(1,nrow(veg_post)), 1)
  
  eta = eta[iter,]
  rho = rho[iter,]
  
} else {
  
  eta = colMeans(eta)
  rho = colMeans(rho)
  
}

eta = unname(eta)[1:K]
rho = unname(rho)[1:K]

# ##########################################################################################################################
# ## chunk: qr decompose X
# ##########################################################################################################################
# 
# x = matrix(1, nrow=(N*T), ncol=1)
# N_p = N*T
# 
# temp = qr(x)
# Q = qr.Q(temp)
# R = qr.R(temp)
# 
# P = Q %*% t(Q)
# # M = diag(N_p) - P
# 
# if (all(P-P[1,1]<1.0e-12)){
#   P = P[1,1]
#   N_p = 1
# }

##########################################################################################################################
## save the data; rdata more efficient, use for processing
##########################################################################################################################
if (kernel == 'gaussian'){ suff = paste0('G_', suff) } else if (kernel == 'pl'){suff = paste0('PL_', suff)}
if (!draw) suff = paste0(suff, '_mean')

dirName = paste0('runs/', N_knots, 'knots_', tmin, 'to', tmax, 'ybp_', suff)

if (one_time){
  dirName = paste0('runs/space_slices_', suff)
}
if (AR){
  dirName = paste0(dirName, '_ar')
}
if (!(file.exists(dirName))) {
  dir.create(dirName)
}


if (one_time){
  subDir=paste0('slice', tmin, 'to', tmax)
  if (!(file.exists(file.path(dirName, subDir)))) {
    dir.create(file.path(dirName, subDir))
  }
} else {
  subDir=paste0('run', dr)
  if (!(file.exists(file.path(dirName, subDir)))) {
    dir.create(file.path(dirName, subDir))
  }
}

fname = file.path(dirName, subDir, 'input')

# note that w is column-major 
save(K, N, T, N_cores, N_knots, res,
     gamma, phi, rho, eta,
     y, 
     idx_cores, 
     d_knots, d_inter, w, #d_pol, #d, 
     lag,
     #        P, N_p, sum_w_pot,
     meta_pol, meta_pol_all,
     sum_w_pot, 
     #pollen_check, # FS - Doesn't exist still...
     knot_coords,
     centers_pls, centers_veg, centers_pol, taxa, ages, y_veg, N_pls,
     file=paste0(fname, '.rdata'))

# convert to row-major # FS - DOESN"T EXIST
if (KW){
  w_new = vector(length=0)
  for (k in 1:K)
    w_new = c(w_new, as.vector(w[k,,]))
  w = array(w_new, c(K, N_cores, N))  
}

dump(c('K', 'N', 'T', 'N_cores', 'N_knots', 'res',
       'gamma', 'phi', 'rho', 'eta',
       'y', 
       'idx_cores', 
       'd_knots', 'd_inter', 'w', #'d_pol', #'d', 
       'lag',
       'sum_w_pot'),
     file=paste0(fname, '.dump'))

##########################################################################################################################
## write meta file with paths
##########################################################################################################################

if (dr==1){
  paths = list(path_grid   = path_grid, 
               path_pls    = path_pls, 
               path_pollen = path_pollen, 
              # path_ages  = path_ages,  # FS - Doesn't exist
               path_cal    = path_cal, 
               path_veg_data = path_veg_data, 
               path_veg_pars = path_veg_pars)
  
  conn=file(file.path(dirName, 'meta.txt'), 'wt')
  write("## Path names", conn)
  for (j in 1:length(paths)) {
    write(paste0('## ', names(paths)[[j]], '=', paths[[j]]), conn)
  }
  close(con=conn)
}
