library(data.table)
library(dplyr)
library(fst)
library(purrr)
library(furrr)
library(stringr)
library(tidyr)
library(rsample)
library(extraDistr)

future::plan("multisession", workers = 6) # Set to desired number of cores


source(file.path("scripts", "99_bayes_funs.R"))



# Settings ----------------------------------------------------------------

# Model priors
mu_0 <- 0.3 # mean (b)
kappa_0 <- 1 # number of observations/precision (c)
a_0 <- 3 # shape of Gamma (g)
b_0 <- 0.2 # rate of Gamma (h)



# Thresholds for inclusion
min_reps <- 3
max_reps <- 25
min_obs <- 30


# Training/test split
train_prop <- 0.8

# Predictions -------------------------------------------------------------

courses <- c("Grandes Lignes", "Stepping Stones")

map(courses, function (course) {
  
  message(course)
  
  # Load data ---------------------------------------------------------------
  
  rof <- read_fst(file.path("data", paste0("rate_of_forgetting_", str_replace_all(course, " ", "_"), ".fst")))
  setDT(rof)
  
  # Only include observations within the specified bounds
  rof_sub <- rof[data.table::between(n_reps, min_reps, max_reps, incbounds = TRUE)]
  
  rm(rof)
  gc()
  
  
  # Training / testing split ------------------------------------------------
  rof_sub_split <- initial_split(rof_sub, prop = train_prop)
  rof_sub_train <- training(rof_sub_split)
  rof_sub_test <- testing(rof_sub_split)
  
  write_fst(rof_sub_train, file.path("data", paste0("training_", str_replace_all(course, " ", "_"), ".fst")))
  write_fst(rof_sub_test, file.path("data", paste0("testing_", str_replace_all(course, " ", "_"), ".fst")))

  message(paste("Total trials:", scales::comma(sum(rof_sub$n_reps))))
  message(paste("Test set trials:", scales::comma(sum(rof_sub_test$n_reps))))
  
  message(paste("Total learning sequences:", scales::comma(nrow(rof_sub))))
  message(paste("Test set learning sequences:", scales::comma(nrow(rof_sub_test))))
  
  
  # Fact-level predictions --------------------------------------------------
  
  message("Fact-level predictions...")
  message(paste("Total facts:", scales::comma(length(unique(rof_sub$fact_id)))))
  message(paste("Test set facts:", scales::comma(length(unique(rof_sub_test$fact_id)))))
  
  # Only include facts with a sufficient number of observations
  facts <- rof_sub_train[, fact_id[.N >= min_obs], by = .(fact_id)][, .(fact_id)]
  
  # Process the data in manageable chunks
  chunk_size <- 1e3
  facts_chunks <- c(seq(1, nrow(facts), by = chunk_size), nrow(facts)+1)
  
  fact_pred <- map_dfr(1:(length(facts_chunks)-1), function (i) {
    
    facts_i <- facts[facts_chunks[i]:(facts_chunks[i+1]-1)]
    
    rof_sub_train[facts_i, on = .(fact_id)] %>%
      nest(data = c(user_id, n_reps, alpha)) %>%
      mutate(res = map(data, ~ run_bayes_model(.$alpha))) %>%
      select(-data) %>%
      unnest(res) %>%
      select(fact_id,
             mu = mu_n,
             kappa = kappa_n,
             a = a_n,
             b = b_n)
    
  })
  
  setDT(fact_pred)
  
  write_fst(fact_pred, file.path("data", "predictions", paste0("pred_fact_", str_replace_all(course, " ", "_"), ".fst")))

  # Compare predicted vs observed values in the test set
  fact_pred_v_obs <- rof_sub_test[fact_pred[, .(fact_id, pred_fact = mu, n_train_obs = kappa - kappa_0)], on = "fact_id"]
  write_fst(fact_pred_v_obs, file.path("data", "predictions", paste0("pred_v_obs_fact_", str_replace_all(course, " ", "_"), ".fst")))
  
  message("Done!")
  
  
  
  # Learner-level predictions -----------------------------------------------
  
  message("Learner-level predictions...")
  message(paste("Total learners:", scales::comma(length(unique(rof_sub$user_id)))))
  message(paste("Test set learners:", scales::comma(length(unique(rof_sub_test$user_id)))))
  
  
  # Only include learners with a sufficient number of observations
  users <- rof_sub_train[, user_id[.N >= min_obs], by = .(user_id)][, .(user_id)]
  
  # Process the data in manageable chunks
  chunk_size <- 1e3
  users_chunks <- c(seq(1, nrow(users), by = chunk_size), nrow(users)+1)
  
  user_pred <- map_dfr(1:(length(users_chunks)-1), function (i) {
    
    users_i <- users[users_chunks[i]:(users_chunks[i+1]-1)]
    
    rof_sub_train[users_i, on = .(user_id)] %>%
      nest(data = c(fact_id, n_reps, alpha)) %>%
      mutate(res = map(data, ~ run_bayes_model(.$alpha))) %>%
      select(-data) %>%
      unnest(res) %>%
      select(user_id,
             mu = mu_n,
             kappa = kappa_n,
             a = a_n,
             b = b_n)
    
  })
  
  setDT(user_pred)
  
  write_fst(user_pred, file.path("data", "predictions", paste0("pred_user_", str_replace_all(course, " ", "_"), ".fst")))

  # Compare predicted vs observed values in the test set
  user_pred_v_obs <- rof_sub_test[user_pred[, .(user_id, pred_user = mu, n_train_obs = kappa - kappa_0)], on = "user_id"]
  write_fst(user_pred_v_obs, file.path("data", "predictions", paste0("pred_v_obs_user_", str_replace_all(course, " ", "_"), ".fst")))
  
  message("Done!")
  
  
  # Fact-and-learner-level predictions --------------------------------------
  
  message("Fact-and-learner-level predictions...")
  
  # Process the data in manageable chunks
  chunk_size <- 1e4
  rof_chunks <- c(seq(1, nrow(rof_sub_test), by = chunk_size), nrow(rof_sub_test)+1)
  
  fact_and_user_pred <- map_dfr(1:(length(rof_chunks)-1), function(i) {
    
    if(i %% 10 == 0) message(scales::percent(i / (length(rof_chunks)-1)))
    
    rof_i <- rof_sub_test[rof_chunks[i]:(rof_chunks[i+1]-1)]
    rof_lst <- split(rof_i, by = c("user_id", "fact_id"), drop = TRUE)
    
    future_map_dfr(rof_lst, function (x) {
      
      f_i <- fact_pred[fact_pred$fact_id == x$fact_id,]
      u_i <- user_pred[user_pred$user_id == x$user_id,]
      
      if(nrow(f_i) == 0 || nrow(u_i) == 0) {
        fact_and_user_pred <- NA
      } else {
        
        fact_t_distr <- calculate_t_distr(f_i$mu, f_i$kappa, f_i$a, f_i$b)
        user_t_distr <- calculate_t_distr(u_i$mu, u_i$kappa, u_i$a, u_i$b)
        
        fact_and_user_pool <- calculate_logarithmic_pool(fact_t_distr, user_t_distr)
        fact_and_user_pred <- get_mode(fact_and_user_pool)
      }
      
      list(fact_id = x$fact_id,
           user_id = x$user_id,
           alpha = x$alpha,
           pred_fact_user = fact_and_user_pred
      )
      
    }, .progress = interactive())
    
  })
  
  setDT(fact_and_user_pred)
  
  write_fst(fact_and_user_pred, file.path("data", "predictions", paste0("pred_fact_and_user_", str_replace_all(course, " ", "_"), ".fst")))
  
  message("Done!")
  
  
})