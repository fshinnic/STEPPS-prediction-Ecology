# from stepps-prediction
# imported 1/6/2026
# to be read by r/pred_process_v1.R
# established functions to read stan binary outputs (I think)

read_stanbin0 <- function(fname, burn=-1) {
  bin = file(fname, "rb")
  nwarmup  = readBin(bin, "integer")
  nsamples = readBin(bin, "integer")
  nparams  = readBin(bin, "integer")
  warmup   = readBin(bin, "numeric", n=nwarmup*nparams, size=4)
  samples  = readBin(bin, "numeric", n=nsamples*nparams, size=4)
  params   = readBin(bin, "character", n=nparams)
  short    = readBin(bin, "character", n=nparams)
  close(bin)
  
  dim(warmup) <- c(nparams, nwarmup)
  dim(samples) <- c(nparams, nsamples)
  
  rownames(warmup) <- params
  rownames(samples) <- params
  
  if (burn > -1 && burn < nsamples) {
    warmup  = NULL
    samples = samples[,(burn+1):nsamples]
  } else {
    warmup = t(warmup)
  }
  
  list(params=params, warmup=warmup, samples=t(samples), par_names=short, nsamples=nsamples, stepsize=0, diagonal=NULL)
}

read_stanbin1 <- function(fname, burn=-1) {
  bin = file(fname, "rb")
  magic    = readBin(bin, "integer")  
  nwarmup  = readBin(bin, "integer")
  nsamples = readBin(bin, "integer")
  nparams  = readBin(bin, "integer")
  ndiag    = readBin(bin, "integer")
  warmup   = readBin(bin, "numeric", n=nwarmup*nparams, size=4)
  stepsize = readBin(bin, "numeric", n=1, size=4)
  diagonal = readBin(bin, "numeric", n=ndiag, size=4)
  samples  = readBin(bin, "numeric", n=nsamples*nparams, size=4)
  params   = readBin(bin, "character", n=nparams)
  short    = readBin(bin, "character", n=nparams)
  close(bin)
  
  dim(warmup) <- c(nparams, nwarmup)
  dim(samples) <- c(nparams, nsamples)
  
  rownames(warmup) <- params
  rownames(samples) <- params
  
  if (burn > -1 && burn < nsamples) {
    warmup  = NULL
    samples = samples[,(burn+1):nsamples]
  } else {
    warmup = t(warmup)
  }
  
  list(params=params, warmup=warmup, samples=t(samples), par_names=short, nsamples=nsamples, stepsize=stepsize, diagonal=diagonal)
}

read_stanbin <- function(fname, burn=-1) {
  bin = file(fname, "rb")
  magic = readBin(bin, "integer")
  if (magic > 0) {
    close(bin)
    read_stanbin0(fname, burn)
  } else if (magic == -1) {
    close(bin)
    read_stanbin1(fname, burn)
  }
}


g_samples_contiguous <- function(dir_name) {
  fname = paste0(dir_name, '/output.bin')
  bin <- read_stanbin(fname)
  
  params = bin$params[which(bin$par_names=="g")]
  samples = t(bin$samples[,which(bin$par_names=="g")])
  ng = nrow(samples)
  
  gname = paste0(dir_name, '/g.bin')
  bin = file(gname, "wb")
  writeBin(nrow(samples), bin) # number of g params
  writeBin(ncol(samples), bin) # number of samples
  for (r in 1:nrow(samples)) {
    writeBin(samples[r,], bin, size=4)
  }
  close(bin)
}