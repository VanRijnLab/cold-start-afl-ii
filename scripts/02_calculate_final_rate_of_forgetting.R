library(data.table)
library(fst)
library(purrr)
library(furrr)
library(stringr)

future::plan("multisession", workers = 6) # Set to desired number of cores


source(file.path("scripts", "99_slimstampen_model_funs.R"))

courses <- c("Grandes Lignes", "Stepping Stones")


map(courses, function (course) {
  
  message(course)
  
  # Load data ---------------------------------------------------------------
  
  d_full <- read_fst(file.path("data", paste0("formatted_", str_replace_all(course, " ", "_"), ".fst")))
  setDT(d_full)
  d_full <- unique(d_full)
  
  # The way we merge facts results in a lot of NA fact IDs -- facts that were never answered correctly
  d <- d_full[!is.na(fact_id)]
  
  # How many trials are discarded because of this?
  message(paste("Total trials:", scales::comma(nrow(d_full))))
  message(paste("Trials discarded because facts were unidentifiable:", scales::comma(nrow(d_full) - nrow(d))))
  
  rm(d_full)
  gc()
  
  setorder(d, user_id, fact_id, start_time)
  
  
  # Calculate rate of forgetting --------------------------------------------
  
  obs <- unique(d[, .(user_id, fact_id)])
  
  # Process the data in manageable chunks
  chunk_size <- 1e4
  obs_chunks <- c(seq(1, nrow(obs), by = chunk_size), nrow(obs)+1)
  
  for(i in 1:(length(obs_chunks)-1)) {
    
    message(paste(i, "/", length(obs_chunks)-1))
    
    obs_i <- obs[obs_chunks[i]:obs_chunks[i+1]-1]
    
    d_i <- d[obs_i, on = .(user_id, fact_id)]
    d_i <- d_i[, .(user_id,
                   fact_id,
                   text = "",
                   start_time,
                   rt,
                   correct,
                   threshold = -0.8)]
    
    d_i_split <- split(d_i, by = c("user_id", "fact_id"), drop = TRUE)
    
    rof_i <- future_map_dfr(d_i_split, function (x) {
      
      n_reps <- nrow(x)
      
      if (n_reps < 3) {
        obs_alpha <- NA
      } else {
        obs_alpha <- calculate_alpha(
          time = max(x$start_time) + 1,
          id = x$fact_id[1],
          factalpha = 0.3,
          responses = x
        )
      }
      
      list(user_id = x$user_id[1],
           fact_id = x$fact_id[1],
           n_reps = n_reps,
           alpha = obs_alpha)
    },
    .progress = interactive()
    )
    
    setDT(rof_i)
    
    write_fst(rof_i, file.path("data", "rof", paste0("rate_of_forgetting_", str_replace_all(course, " ", "_"), "_", i, ".fst")))
    
  }
  
  
  # Combine files -----------------------------------------------------------
  
  rof_files <- list.files(file.path("data", "rof"), 
                          pattern = paste0(".*", str_replace_all(course, " ", "_"), ".*"),
                          full.names = TRUE)
  
  rof <- map_dfr(rof_files, read_fst)
  write_fst(rof, file.path("data", paste0("rate_of_forgetting_", str_replace_all(course, " ", "_"), ".fst")))
  
})