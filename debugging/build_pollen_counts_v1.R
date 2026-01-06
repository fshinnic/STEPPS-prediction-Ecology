build_pollen_counts <- function(tmin, tmax, int,
                                pollen_ts,
                                taxa_all,
                                taxa_sub,
                                age_model) {
  
  ## ---- choose age column ----
  if (age_model == "bacon") {
    age_col <- "age_bacon"
  } else if (age_model == "bchron") {
    age_col <- "age_bchron"
  } else {
    age_col <- "age_default"
  }
  
  ## ---- find first taxa column ----
  taxa.start.col <- min(match(taxa_all, colnames(pollen_ts)), na.rm = TRUE)
  
  ## ---- subset to time window (int) ----
  keep <- pollen_ts[, age_col] >= tmin & pollen_ts[, age_col] <= tmax
  pollen_sub <- pollen_ts[keep, , drop = FALSE]
  
  meta_pol <- pollen_sub[, 1:(taxa.start.col - 1), drop = FALSE]
  counts   <- pollen_sub[, taxa.start.col:ncol(pollen_sub), drop = FALSE]
  
  ## add age placeholder
  meta_pol$age <- NA_real_
  
  ## ---- bin breaks ----
  breaks <- seq(tmin, tmax, by = int)
  
  ## ---- initialize outputs (DATA FRAMES ONLY) ----
  counts_agg <- counts[0, ]
  
  meta_agg <- data.frame(
    meta_pol[0, 1:6],
    age  = numeric(0),
    zero = logical(0)
  )
  
  meta_all <- meta_pol[0, ]
  
  ## ---- loop over cores ----
  ids <- unique(meta_pol$id)
  
  for (i in seq_along(ids)) {
    
    # tells when summing each core
    message("core ", i)
    
    core_rows <- which(meta_pol$id == ids[i])
    
    for (j in seq_len(length(breaks) - 1)) {
      
      age_mid <- breaks[j] + int / 2
      
      age_rows <- core_rows[
        meta_pol[core_rows, age_col] >= breaks[j] &
          meta_pol[core_rows, age_col] <  breaks[j + 1]
      ]
      
      ## ---- cases ----
      if (length(age_rows) > 1) {
        
        counts_row <- as.data.frame(
          t(colSums(counts[age_rows, , drop = FALSE]))
        )
        
        meta_row <- data.frame(
          meta_pol[age_rows[1], 1:6, drop = FALSE],
          age  = age_mid / 100,
          zero = FALSE
        )
        
        meta_all_rows <- meta_pol[age_rows, , drop = FALSE]
        meta_all_rows$age <- age_mid / 100
        
      } else if (length(age_rows) == 1) {
        
        counts_row <- counts[age_rows, , drop = FALSE]
        
        meta_row <- data.frame(
          meta_pol[age_rows, 1:6, drop = FALSE],
          age  = age_mid / 100,
          zero = FALSE
        )
        
        meta_all_rows <- meta_pol[age_rows, , drop = FALSE]
        meta_all_rows$age <- age_mid / 100
        
      } else {
        
        counts_row <- as.data.frame(
          t(rep(0, ncol(counts)))
        )
        colnames(counts_row) <- colnames(counts)
        
        meta_row <- data.frame(
          meta_pol[core_rows[1], 1:6, drop = FALSE],
          age  = age_mid / 100,
          zero = TRUE
        )
        
        meta_all_rows <- meta_pol[core_rows[1], , drop = FALSE]
        meta_all_rows$age <- NA_real_
      }
      
      ## ---- bind ----
      counts_agg <- rbind(counts_agg, counts_row)
      meta_agg   <- rbind(meta_agg, meta_row)
      meta_all   <- rbind(meta_all, meta_all_rows)
    }
  }
  
  return(list(
    counts_agg = counts_agg,
    meta_agg   = meta_agg,
    meta_all   = meta_all
  ))
}
