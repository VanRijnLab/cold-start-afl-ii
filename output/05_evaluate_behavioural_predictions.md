Evaluate behavioural predictions
================
Maarten van der Velde
Last updated: 2022-10-17

-   [Overview](#overview)
-   [Setup](#setup)
    -   [Helper functions](#helper-functions)
-   [Calculate predictions](#calculate-predictions)
-   [Activation](#activation)
-   [Response time](#response-time)
    -   [Predicted response time](#predicted-response-time)
        -   [Distribution of predictions](#distribution-of-predictions)
        -   [Predicted vs observed values](#predicted-vs-observed-values)
-   [Response accuracy](#response-accuracy)
    -   [Predicted response accuracy](#predicted-response-accuracy)
        -   [Distribution of predictions](#distribution-of-predictions-1)
        -   [Predicted vs observed values](#predicted-vs-observed-values-1)
-   [Combined plot](#combined-plot)
-   [Session info](#session-info)

# Overview

This notebook evaluates the behavioural predictions (RT, accuracy) that follow from the predicted rates of forgetting. We simulate the behavioural predictions that the adaptive fact learning system would have made, if it had used the predicted rate of forgetting as the starting estimate in the learning sequence.

# Setup

``` r
library(fst)
library(data.table)
library(tidyr)
library(purrr)
library(furrr)
library(stringr)
library(ggplot2)
library(patchwork)
library(wesanderson)
library(lme4)
library(lmerTest)
library(multcomp)
```

``` r
source(file.path("..", "scripts", "99_slimstampen_model_funs.R"))
```

``` r
future::plan("multiprocess", workers = 6) # Set to desired number of cores
```

``` r
theme_set(theme_light(base_size = 14) +
            theme(strip.text = element_text(colour = "black")))

condition_colours <- wes_palette("Darjeeling1", n = 5)
condition_colours[c(2, 4, 5)] <- condition_colours[c(4, 5, 2)]

dataset_colours <- wes_palette("Darjeeling2", n = 5)[c(2, 3)]
```

## Helper functions

``` r
load_data_with_predictions <- function (course) {
  
  # Data
  d_full <- read_fst(file.path("..", "data", paste0("formatted_", str_replace_all(course, " ", "_"), ".fst")))
  setDT(d_full)
  d <- d_full[!is.na(fact_id), .(user_id, fact_id, start_time, rt, correct)]
  rm(d_full)
  gc()
  setorder(d, user_id, fact_id, start_time)
  
  # ROF predictions
  pred_user <- read_fst(file.path("..", "data", "predictions", paste0("pred_v_obs_user_", str_replace_all(course, " ", "_"), ".fst")))
  pred_fact <- read_fst(file.path("..", "data", "predictions", paste0("pred_v_obs_fact_", str_replace_all(course, " ", "_"), ".fst")))
  pred_fact_user <- read_fst(file.path("..", "data", "predictions", paste0("pred_fact_and_user_", str_replace_all(course, " ", "_"), ".fst")))
  setDT(pred_user)
  setDT(pred_fact)
  setDT(pred_fact_user)
  
  # Remove NA and duplicates
  pred_user <- unique(pred_user[!is.na(alpha)])
  pred_fact <- unique(pred_fact[!is.na(alpha)])
  pred_fact_user <- unique(pred_fact_user[!is.na(alpha)])
  
  # Make Domain prediction
  pred_domain <- mean(unique(pred_fact, by = c("fact_id"))$pred_fact)
  pred_default <- 0.3
  
  # Combine
  setnames(pred_user, "n_train_obs", "n_train_obs_user")
  setnames(pred_fact, "n_train_obs", "n_train_obs_fact")
  pred_all <- merge(pred_user, pred_fact, by = c("user_id", "fact_id", "alpha", "n_reps"), all = TRUE)
  pred_all <- merge(pred_all, pred_fact_user, by = c("user_id", "fact_id", "alpha"), all = TRUE)
  pred_all[, pred_default := pred_default]
  pred_all[, pred_domain := pred_domain]

  d_prep <- d[, .(user_id,
                  fact_id,
                  text = "",
                  start_time,
                  rt,
                  correct,
                  threshold = -0.8)]
  
  d_pred <- merge(pred_all, d_prep, by = c("user_id", "fact_id"))
  
  return(d_pred)
}

predict_behaviour <- function (d) {
  
    # Process the data in manageable chunks
    chunk_size <- 1e4
    obs <- unique(d[, .(user_id, fact_id)])
    obs_chunks <- c(seq(1, nrow(obs), by = chunk_size), nrow(obs)+1)
    
    pred_beh <- map_dfr(1:(length(obs_chunks)-1), function (i) {
      
      msg <- paste("Chunk", i, "/", length(obs_chunks)-1)
      system(paste("echo", msg))
      
      obs_i <- obs[obs_chunks[i]:obs_chunks[i+1]-1]
      d_i <- d[obs_i, on = .(user_id, fact_id)]
      d_i_list <- split(d_i, by = c("user_id", "fact_id"), drop = TRUE)
      
      pred_beh_i <- future_map_dfr(d_i_list, function (learn_seq) {
        
        # Organise predictions for this sequence
        learn_seq_preds <- pivot_longer(
          learn_seq[1,], 
          pred_user:pred_domain,
          names_to = "prediction_type",
          names_prefix = "pred_",
          values_to = "predicted_alpha"
        )
        setDT(learn_seq_preds)
        
        # Look only at trial 3 in each sequence
        trial <- learn_seq[3,]
        delay <- learn_seq[2:3, diff(start_time)]
        
        # Calculate behavioural predictions for each prediction method
        map_dfr(seq(nrow(learn_seq_preds)), function (j) {
          
          prediction_type <- learn_seq_preds[j, prediction_type]
          predicted_alpha <- learn_seq_preds[j, predicted_alpha]
    
          predicted_activation <- NA
          predicted_accuracy <- NA
          predicted_rt <- NA
    
          if (!is.na(predicted_alpha)) {
            
            predicted_activation <- calculate_activation(
              time = trial[, start_time],
              id = trial[, fact_id],
              factalpha = predicted_alpha,
              responses = learn_seq[1:2,]
            )
      
            predicted_accuracy <- p_recall(
              activation = predicted_activation,
              threshold = -0.8,
              activation_noise = 0.5
            )
            
            predicted_rt <- estimate_reaction_time_from_activation(
              activation = predicted_activation,
              reading_time = 300
            )
          }
          
          return(
            list(
              user_id = trial[, user_id],
              fact_id = trial[, fact_id],
              delay = delay,
              correct = trial[, correct],
              rt = trial[, rt],
              prediction_type = prediction_type,
              predicted_alpha = predicted_alpha,
              predicted_activation = predicted_activation,
              predicted_accuracy = predicted_accuracy,
              predicted_rt = predicted_rt
            )
          )
          
        })
      })
    })
    
    setDT(pred_beh)
    
    # Remove rows without prediction
    pred_beh <- pred_beh[!is.na(predicted_alpha)]
    pred_beh <- pred_beh[!is.infinite(predicted_rt)]
    
    # Set proper condition labels
    condition_labels <- data.table(
      prediction_type = c("default", "domain", "fact", "user", "fact_user"),
      prediction_label = factor(
        c("Default", "Domain", "Fact", "Learner", "Fact & Learner"),
        levels = c("Default", "Domain", "Fact", "Learner", "Fact & Learner")
      )
    )
    pred_beh <- pred_beh[condition_labels, on = .(prediction_type)]
    
    
    return(pred_beh)
}
```

# Calculate predictions

Load test set data with rate of forgetting predictions:

``` r
pred_gl <- load_data_with_predictions("Grandes Lignes")
pred_ss <- load_data_with_predictions("Stepping Stones")
```

Calculate behavioural predictions for trial 3:

``` r
pred_gl_beh_path <- file.path("..", "data", "predictions", "pred_behaviour_gl.fst")
pred_ss_beh_path <- file.path("..", "data", "predictions", "pred_behaviour_ss.fst")

if (!file.exists(pred_gl_beh_path)) {
  pred_gl_beh <- predict_behaviour(pred_gl)
  write_fst(pred_gl_beh, pred_gl_beh_path)
} else {
  pred_gl_beh <- read_fst(pred_gl_beh_path)
  setDT(pred_gl_beh)
}

if (!file.exists(pred_ss_beh_path)) {
  pred_ss_beh <- predict_behaviour(pred_ss)
  write_fst(pred_ss_beh, pred_ss_beh_path)
} else {
  pred_ss_beh <- read_fst(pred_ss_beh_path)
  setDT(pred_ss_beh)
}
```

``` r
pred_beh <- rbind(pred_gl_beh[, course := "French"],
                  pred_ss_beh[, course := "English"])

rm(pred_gl_beh, pred_ss_beh)
gc()
```

    ##              used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells    3651434  195.1    7217305   385.5    7217305   385.5
    ## Vcells 1237212008 9439.2 2778101033 21195.3 2892276405 22066.4

# Activation

``` r
p_act_dist <- ggplot(pred_beh[between(predicted_activation, -2, 0)], aes(x = predicted_activation, fill = prediction_label)) +
  facet_grid(course ~ prediction_label, scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .01) +
  guides(fill = "none") +
  labs(x = "RT (in s)",
       y = "Density") +
  scale_fill_manual(values = condition_colours)

# p_act_dist

ggsave(plot = p_act_dist, file.path("..", "output", "activation_distribution.png"),
       device = "png", width = 7.5, height = 4.5)

rm(p_act_dist)
```

# Response time

Distribution of observed correct RT (truncated at 25 seconds for readability):

``` r
p_rt_dist <- ggplot(pred_beh[correct == 1 & between(rt, 0, 25000)], aes(x = rt/1000)) +
  facet_grid(course ~ ., scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .1) +
  guides(fill = "none") +
  labs(x = "RT (in s)",
       y = "Density")

# p_rt_dist

ggsave(plot = p_rt_dist, file.path("..", "output", "rt_observed_distribution.png"),
       device = "png", width = 6, height = 6)

rm(p_rt_dist)
```

## Predicted response time

### Distribution of predictions

``` r
p_rt_pred_dist <- ggplot(pred_beh, aes(x = predicted_rt/1000)) +
  facet_grid(course ~ ., scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .1) +
  guides(fill = "none") +
  labs(x = "Predicted RT (in s)",
       y = "Density")

# p_rt_pred_dist

ggsave(plot = p_rt_pred_dist, file.path("..", "output", "rt_predicted_distribution.png"),
       device = "png", width = 6, height = 6)

rm(p_rt_pred_dist)
```

### Predicted vs observed values

Calculate absolute prediction error:

``` r
pred_rt_error <- (pred_beh
                  [correct == 1]
                  [between(rt, 0, 25000)]
                  [, rt_pred_error := predicted_rt - rt]
                  [, abs_rt_pred_error := abs(rt_pred_error)]
)

pred_rt_error_avg <- pred_rt_error[, .(mae = mean(abs_rt_pred_error),
                                       ae_se = sd(abs_rt_pred_error)/.N), 
                                   by = .(course, prediction_label)]
  
n_obs <- pred_rt_error[, .N, by = .(course, prediction_label)]
```

``` r
# plot_range <- range(pred_beh$predicted_rt/1000, na.rm = TRUE)
plot_range <- c(0, 10)
plot_breaks <- seq(0, 10, by = 2)
# 
# plot_range <- quantile(pred_beh$predicted_rt/1000, c(.005, .995))

p_rt_pred_v_obs <- ggplot(pred_beh[correct == 1], aes(x = predicted_rt/1000, y = rt/1000, colour = prediction_label)) +
  facet_grid(course ~ prediction_label) +
  geom_abline(slope = 1, intercept = 0, lty = 3, alpha = 0.75) +
  geom_point(alpha = .1, size = .1, pch = ".") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
  geom_label(data = pred_rt_error_avg,
            aes(label = paste("MAE =", formatC(mae/1000, digits = 3, flag = "#"))),
            x = plot_range[2], y = plot_range[1],
            hjust = 1, colour = "NA", size = 3,
            alpha = .9,
            label.size = NA) +
  geom_text(data = pred_rt_error_avg,
            aes(label = paste("MAE =", formatC(mae/1000, digits = 3, flag = "#"))),
            x = plot_range[2], y = plot_range[1],
            hjust = 1, colour = "black", size = 3) +
  geom_label(data = n_obs,
            aes(label = paste("n =", scales::comma(N))),
            x = plot_range[2],
            y = plot_range[2],
            hjust = 1, colour = "NA", size = 3,
            alpha = .9,
            label.size = NA) +
  geom_text(data = n_obs,
            aes(label = paste("n =", scales::comma(N))),
            x = plot_range[2],
            y = plot_range[2],
            hjust = 1, colour = "black", size = 3) +
  guides(colour = "none") +
  labs(x = "Predicted RT (s)",
       y = "Observed RT (s)") +
  coord_fixed(ratio = 1, xlim = plot_range, ylim = plot_range) +
  scale_x_continuous(breaks = plot_breaks) +
  scale_y_continuous(breaks = plot_breaks) +
  scale_colour_manual(values = condition_colours)

p_rt_pred_v_obs
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
ggsave(plot = p_rt_pred_v_obs, file.path("..", "output", "rt_predicted_vs_observed.png"),
       device = "png", width = 10, height = 4.5)

rm(p_rt_pred_v_obs)
```

#### Prediction error

Distribution of prediction error (truncated to \[-5, 5\] for readability):

``` r
p_rt_pred_error <- ggplot(pred_rt_error, aes(x = rt_pred_error/1000, fill = prediction_label)) +
  facet_grid(prediction_label ~ course , scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .1) +
  guides(fill = "none") +
  labs(x = "RT prediction error in s (predicted - observed)",
       y = "Density") +
  coord_cartesian(xlim = c(-5, 5)) +
  scale_fill_manual(values = condition_colours)

p_rt_pred_error
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
ggsave(plot = p_rt_pred_error, file.path("..", "output", "rt_prediction_error.png"),
       device = "png", width = 5, height = 7.5)

rm(p_rt_pred_error)
```

#### Absolute prediction error

``` r
p_rt_abs_pred_error <- ggplot(pred_rt_error, aes(x = abs_rt_pred_error/1000, fill = prediction_label)) +
  facet_grid(prediction_label ~ course, scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .1) +
  guides(fill = "none") +
  labs(x = "RT prediction error in s (predicted - observed)",
       y = "Density") +
  coord_cartesian(xlim = c(0, 5)) +
  scale_fill_manual(values = condition_colours)

p_rt_abs_pred_error
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
ggsave(plot = p_rt_abs_pred_error, file.path("..", "output", "rt_absolute_prediction_error.png"),
       device = "png", width = 5, height = 7.5)

rm(p_rt_abs_pred_error)
```

``` r
ggplot(pred_rt_error_avg, aes(x = prediction_label, y = mae/1000, colour = course)) +
  geom_boxplot(data = pred_rt_error,
               aes(y = abs_rt_pred_error/1000, group = interaction(course, prediction_label)),
               colour = "grey70",
               width = .25,
               outlier.shape = NA,
               position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = mae/1000 - ae_se/1000, ymax = mae/1000 + ae_se/1000), width = 0, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(x = "Method",
       y = "Absolute RT prediction error (in s)",
       colour = "Course")
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-16-1.png)

Fit a regression model on absolute RT prediction error.

##### French

``` r
m_rt_pred_error_gl_file <- file.path("..", "data", "model_fits", "m_rt_pred_error_Grandes_Lignes.rda")

if (file.exists(m_rt_pred_error_gl_file)) {
  load(m_rt_pred_error_gl_file)
} else {
  
  pred_gl_reg <- (
    pred_rt_error
    [course == "French"]
    [sample(.N, 1e6)]
    [, .(prediction_label, abs_rt_pred_error, user_id, fact_id)]
  )
  
  m_rt_pred_error_gl <- lmer(abs_rt_pred_error ~ prediction_label + 
                               (1 | user_id) + (1 | fact_id),
                             data = pred_gl_reg,
                             control = lmerControl(optimizer ="bobyqa"))
  
  save(m_rt_pred_error_gl, file = m_rt_pred_error_gl_file)
}

summary(m_rt_pred_error_gl)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## abs_rt_pred_error ~ prediction_label + (1 | user_id) + (1 | fact_id)
    ##    Data: pred_gl_reg
    ## Control: lmerControl(optimizer = "bobyqa")
    ## 
    ## REML criterion at convergence: 18047298
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5110 -0.4053 -0.1761  0.0625 11.1385 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  user_id  (Intercept)  220053   469.1  
    ##  fact_id  (Intercept)  194518   441.0  
    ##  Residual             3843235  1960.4  
    ## Number of obs: 1000000, groups:  user_id, 40820; fact_id, 22762
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      1433.560      6.509  53722.081 220.257
    ## prediction_labelDomain            -19.468      6.192 974447.696  -3.144
    ## prediction_labelFact             -118.481      6.218 974827.243 -19.056
    ## prediction_labelLearner           -42.334      6.270 975733.404  -6.751
    ## prediction_labelFact & Learner   -101.382      6.300 976077.378 -16.091
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## prediction_labelDomain          0.00167 ** 
    ## prediction_labelFact            < 2e-16 ***
    ## prediction_labelLearner        1.46e-11 ***
    ## prediction_labelFact & Learner  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) prdc_D prdc_F prdc_L
    ## prdctn_lblD -0.474                     
    ## prdctn_lblF -0.468  0.497              
    ## prdctn_lblL -0.462  0.492  0.490       
    ## prdctn_lF&L -0.455  0.490  0.489  0.486

Compare different prediction types to each other:

``` r
ht_rt_gl <- glht(m_rt_pred_error_gl, linfct = mcp(prediction_label = "Tukey"))
summary(ht_rt_gl)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = abs_rt_pred_error ~ prediction_label + (1 | user_id) + 
    ##     (1 | fact_id), data = pred_gl_reg, control = lmerControl(optimizer = "bobyqa"))
    ## 
    ## Linear Hypotheses:
    ##                               Estimate Std. Error z value Pr(>|z|)    
    ## Domain - Default == 0          -19.468      6.192  -3.144  0.01418 *  
    ## Fact - Default == 0           -118.481      6.218 -19.056  < 0.001 ***
    ## Learner - Default == 0         -42.334      6.270  -6.751  < 0.001 ***
    ## Fact & Learner - Default == 0 -101.382      6.300 -16.091  < 0.001 ***
    ## Fact - Domain == 0             -99.013      6.225 -15.905  < 0.001 ***
    ## Learner - Domain == 0          -22.867      6.278  -3.642  0.00256 ** 
    ## Fact & Learner - Domain == 0   -81.914      6.308 -12.986  < 0.001 ***
    ## Learner - Fact == 0             76.146      6.304  12.079  < 0.001 ***
    ## Fact & Learner - Fact == 0      17.099      6.329   2.702  0.05375 .  
    ## Fact & Learner - Learner == 0  -59.047      6.374  -9.264  < 0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Inspect the model's residuals:

``` r
qqnorm(resid(m_rt_pred_error_gl))
qqline(resid(m_rt_pred_error_gl), col = "red")
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
plot(m_rt_pred_error_gl)
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-20-1.png)

##### English

``` r
m_rt_pred_error_ss_file <- file.path("..", "data", "model_fits", "m_rt_pred_error_Stepping_Stones.rda")

if (file.exists(m_rt_pred_error_ss_file)) {
  load(m_rt_pred_error_ss_file)
} else {
  
  pred_ss_reg <- (
    pred_rt_error
    [course == "English"]
    [sample(.N, 1e6)]
    [, .(prediction_label, abs_rt_pred_error, user_id, fact_id)]
  )
  
  m_rt_pred_error_ss <- lmer(abs_rt_pred_error ~ prediction_label + 
                               (1 | user_id) + (1 | fact_id),
                             data = pred_ss_reg,
                             control = lmerControl(optimizer ="bobyqa")
  )
  
  save(m_rt_pred_error_ss, file = m_rt_pred_error_ss_file)
}

summary(m_rt_pred_error_ss)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## abs_rt_pred_error ~ prediction_label + (1 | user_id) + (1 | fact_id)
    ##    Data: pred_ss_reg
    ## Control: lmerControl(optimizer = "bobyqa")
    ## 
    ## REML criterion at convergence: 18097553
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6170 -0.3973 -0.1997  0.0215 10.9768 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  user_id  (Intercept)  150645   388.1  
    ##  fact_id  (Intercept)   72213   268.7  
    ##  Residual             4075141  2018.7  
    ## Number of obs: 1000000, groups:  user_id, 85884; fact_id, 45600
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                      1312.875      5.205 161344.300 252.225
    ## prediction_labelDomain            -27.214      6.416 984927.056  -4.242
    ## prediction_labelFact              -58.734      6.448 985192.856  -9.109
    ## prediction_labelLearner           -41.073      6.467 984764.807  -6.351
    ## prediction_labelFact & Learner    -64.932      6.493 985035.484 -10.000
    ##                                Pr(>|t|)    
    ## (Intercept)                     < 2e-16 ***
    ## prediction_labelDomain         2.22e-05 ***
    ## prediction_labelFact            < 2e-16 ***
    ## prediction_labelLearner        2.14e-10 ***
    ## prediction_labelFact & Learner  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) prdc_D prdc_F prdc_L
    ## prdctn_lblD -0.617                     
    ## prdctn_lblF -0.613  0.498              
    ## prdctn_lblL -0.609  0.496  0.494       
    ## prdctn_lF&L -0.606  0.494  0.492  0.491

Compare different prediction types to each other:

``` r
ht_rt_ss <- glht(m_rt_pred_error_ss, linfct = mcp(prediction_label = "Tukey"))
summary(ht_rt_ss)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = abs_rt_pred_error ~ prediction_label + (1 | user_id) + 
    ##     (1 | fact_id), data = pred_ss_reg, control = lmerControl(optimizer = "bobyqa"))
    ## 
    ## Linear Hypotheses:
    ##                               Estimate Std. Error z value Pr(>|z|)    
    ## Domain - Default == 0          -27.214      6.416  -4.242   <0.001 ***
    ## Fact - Default == 0            -58.734      6.448  -9.109   <0.001 ***
    ## Learner - Default == 0         -41.073      6.467  -6.351   <0.001 ***
    ## Fact & Learner - Default == 0  -64.932      6.493 -10.000   <0.001 ***
    ## Fact - Domain == 0             -31.520      6.445  -4.890   <0.001 ***
    ## Learner - Domain == 0          -13.859      6.466  -2.143   0.2017    
    ## Fact & Learner - Domain == 0   -37.717      6.491  -5.811   <0.001 ***
    ## Learner - Fact == 0             17.661      6.497   2.718   0.0514 .  
    ## Fact & Learner - Fact == 0      -6.198      6.521  -0.950   0.8770    
    ## Fact & Learner - Learner == 0  -23.859      6.540  -3.648   0.0025 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Inspect the model's residuals:

``` r
qqnorm(resid(m_rt_pred_error_ss))
qqline(resid(m_rt_pred_error_ss), col = "red")
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
plot(m_rt_pred_error_ss)
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-24-1.png)

##### Comparison

``` r
ht_rt_gl_tidy <- broom::tidy(confint(ht_rt_gl))
ht_rt_ss_tidy <- broom::tidy(confint(ht_rt_ss))
setDT(ht_rt_gl_tidy)
setDT(ht_rt_ss_tidy)

ht_rt_both_tidy <- rbind(ht_rt_gl_tidy[, course := "French"],
                      ht_rt_ss_tidy[, course := "English"])
```

``` r
p_rt_pred_error_comp <- ggplot(ht_rt_both_tidy, aes(x = lhs, y = estimate, ymin = conf.low, ymax = conf.high, colour = course)) +
  geom_hline(yintercept = 0, linetype = "11", colour = "grey60") +
  geom_errorbar(width = 0.1) + 
  geom_point() +
  labs(x = "Linear hypotheses",
       y = "Estimate",
       caption = "Tukey's range test. Error bars show 95% family-wise confidence level.",
       colour = "Course") +
  coord_flip()

p_rt_pred_error_comp
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
ggsave(plot = p_rt_pred_error_comp, file.path("..", "output", "rt_prediction_error_comparisons.png"),
       device = "png", width = 7.5, height = 5)

rm(p_rt_pred_error_comp)
```

##### Summary plot

``` r
pred_rt_error_avg[, prediction_rank := frank(-mae), by = .(course)]

annotation_df_ss <- data.table(
  course = rep("English", 10),
  start = c(1, 1, 1, 1,
            2, 2, 2,
            3, 3,
            4
  ),
  end = c(2, 3, 4, 5,
          3, 4, 5,
          4, 5,
          5
  ),
  y = seq(max(pred_rt_error_avg$mae)*1.01 + 45, max(pred_rt_error_avg$mae)*1.01, by = -5),
  label = c("p < .001", "p < .001", "p < .001", "p < .001",
            "n.s.", "p < .001", "p < .001",
            "n.s.", "p < .01",
            "n.s.")
)

annotation_df_gl <- data.table(
  course = rep("French", 10),
  start = c(1, 1, 1, 1,
            2, 2, 2,
            3, 3,
            4
  ),
  end = c(2, 3, 4, 5,
          3, 4, 5,
          4, 5,
          5
  ),
  y = seq(max(pred_rt_error_avg$mae)*1.01 + 45, max(pred_rt_error_avg$mae)*1.01, by = -5),
  label = c("p < .05", "p < .001", "p < .001", "p < .001",
            "p < .01", "p < .001", "p < .001",
            "p < .001", "p < .001",
            "n.s.")
)

annotation_df_rt <- rbind(annotation_df_ss, annotation_df_gl)
annotation_df_rt[, label := factor(label, levels = c("p < .001", "p < .01", "p < .05", "n.s."))]

p_rt_pred_error_summary <- ggplot(pred_rt_error_avg, aes(x = prediction_rank, y = mae)) +
  facet_grid(~ course) +
  geom_line(data = annotation_df_rt,
            aes(x = 1, y = 1300, lty = label, alpha = label, colour = NULL)) + # Dummy line to get legend
  geom_line(aes(colour = course, group = course)) +
  geom_errorbar(aes(ymin = mae - ae_se, ymax = mae + ae_se), width = 0) +
  geom_point(aes(colour = course, group = course)) +
  geom_label(aes(label = prediction_label), 
             colour = "black", 
             alpha = .9,
             label.size = NA, 
             nudge_y = -15) +
  labs(x = NULL,
       y = "Absolute prediction error:\nresponse time (s)",
       colour = "Course") +
  scale_x_continuous(expand = expansion(add = .75), breaks = NULL) +
  scale_y_continuous(labels = scales::comma_format(big.mark = ".")) +
  scale_colour_manual(values = dataset_colours) +
  scale_linetype_manual(values = c("p < .001" = 1,
                                   "p < .01" = 5,
                                   "p < .05" = 2,
                                   "n.s." = 3),
                        name = "Pairwise comparison:") +
  scale_alpha_manual(values = c("p < .001" = 1,
                                "p < .01" = .75,
                                "p < .05" = .5, 
                                "n.s." = .25), 
                     name = "Pairwise comparison:") +
  guides(colour = "none") +
  ggsignif::geom_signif(data = annotation_df_rt,
                        aes(xmin = start, xmax = end, annotations = "", 
                            y_position = y, lty = label, alpha = label),
                        tip_length = 0,
                        manual = TRUE)  +
  theme(legend.position = "bottom",
        legend.justification = "right")
```

    ## Warning: Ignoring unknown aesthetics: xmin, xmax, annotations, y_position

``` r
p_rt_pred_error_summary
```

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
ggsave(file.path("..", "output", "rt_absolute_prediction_error_summary.png"),
       device = "png", width = 10, height = 4)
```

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

##### Improvement

How big was the improvement from worst to best prediction method?

French:

``` r
# Absolute change
ht_rt_gl_tidy[lhs == "Fact - Default", estimate[1]]
```

    ## [1] -118.4807

``` r
# % change
scales::percent(
  ht_rt_gl_tidy[lhs == "Fact - Default", estimate[1]] / fixef(m_rt_pred_error_gl)[[1]],
  accuracy = .1)
```

    ## [1] "-8.3%"

English:

``` r
# Absolute change
ht_rt_ss_tidy[lhs == "Fact & Learner - Default", estimate[1]]
```

    ## [1] -64.93171

``` r
# % change
scales::percent(
  ht_rt_ss_tidy[lhs == "Fact & Learner - Default", estimate[1]] / fixef(m_rt_pred_error_ss)[[1]],
  accuracy = .1)
```

    ## [1] "-4.9%"

# Response accuracy

## Predicted response accuracy

### Distribution of predictions

``` r
p_acc_pred_dist <- ggplot(pred_beh, aes(x = predicted_accuracy, fill = prediction_label)) +
    facet_grid(course ~ prediction_label, scales = "free_y") +
    geom_histogram(binwidth = .01) +
    guides(fill = "none") +
    labs(x = "Predicted accuracy") +
    scale_x_continuous(breaks = seq(0, 1, by = .25), labels = scales::percent_format()) +
    scale_fill_manual(values = condition_colours)
  
p_acc_pred_dist
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
ggsave(plot = p_acc_pred_dist, file.path("..", "output", "acc_predicted_distribution.png"),
         device = "png", width = 6, height = 7.5)

rm(p_acc_pred_dist)
```

### Predicted vs observed values

``` r
plot_dodge <- function(y, dodge = .1) {
  return (y * (1 + dodge) - dodge/2)
}
```

``` r
p_acc_pred_v_obs <- ggplot(pred_beh, aes(x = predicted_accuracy, y = correct, group = prediction_label, colour = prediction_label, fill = prediction_label)) +
    facet_grid(course ~ prediction_label) +
    geom_point(aes(y = correct),
               position = position_jitter(width = 0, height = .025, seed = 123),
               size = .001, pch = ".", alpha = .1) +
    labs(x = "Predicted accuracy",
         y = "Response accuracy",
         colour = "Prediction method",
         fill = "Prediction method") +
  guides(colour = "none",
         fill = "none") +
    scale_x_continuous(breaks = seq(0, 1, by = .25), labels = scales::percent_format()) +
    scale_y_continuous(breaks = seq(0, 1, by = .25), labels = scales::percent_format()) +
    scale_colour_manual(values = condition_colours) +
    scale_fill_manual(values = condition_colours) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")

p_acc_pred_v_obs
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-32-1.png)

``` r
ggsave(plot = p_acc_pred_v_obs, file.path("..", "output", "acc_predicted_vs_observed.png"),
       device = "png", width = 10, height = 4.5)

rm(p_acc_pred_v_obs)
```

#### Prediction error

``` r
pred_acc_error <- (pred_beh
                   [, acc_pred_error := predicted_accuracy - correct]
                   [, abs_acc_pred_error := abs(acc_pred_error)])


pred_acc_error_avg <- pred_acc_error[, .(mae = mean(abs_acc_pred_error),
                                         ae_se = sd(abs_acc_pred_error)/.N), 
                                     by = .(course, prediction_label)]

n_obs <- pred_acc_error[, .N, by = .(course, prediction_label)]
```

Distribution of prediction error:

``` r
p_acc_pred_error <- ggplot(pred_acc_error, aes(x = acc_pred_error, fill = prediction_label)) +
  facet_grid(prediction_label ~ course , scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .01) +
  guides(fill = "none") +
  labs(x = "Accuracy prediction error (predicted - observed)",
       y = "Density") +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_fill_manual(values = condition_colours)

p_acc_pred_error
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-34-1.png)

``` r
ggsave(plot = p_acc_pred_error, file.path("..", "output", "acc_prediction_error.png"),
       device = "png", width = 5, height = 7.5)

rm(p_acc_pred_error)
```

#### Absolute prediction error

``` r
p_abs_acc_pred_error <- ggplot(pred_acc_error, aes(x = abs_acc_pred_error, fill = prediction_label)) +
  facet_grid(prediction_label ~ course , scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .01) +
  guides(fill = "none") +
  labs(x = "Accuracy prediction error (predicted - observed)",
       y = "Density") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_fill_manual(values = condition_colours)

p_abs_acc_pred_error
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-35-1.png)

``` r
ggsave(plot = p_abs_acc_pred_error, file.path("..", "output", "acc_absolute_prediction_error.png"),
       device = "png", width = 5, height = 7.5)

rm(p_abs_acc_pred_error)
```

``` r
ggplot(pred_acc_error_avg, aes(x = prediction_label, y = mae, colour = course)) +
  geom_boxplot(data = pred_acc_error,
               aes(y = abs_acc_pred_error, group = interaction(course, prediction_label)),
               colour = "grey70",
               width = .25,
               outlier.shape = NA,
               position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = mae - ae_se, ymax = mae + ae_se), width = 0, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Method",
       y = "Absolute accuracy prediction error",
       colour = "Course")
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-36-1.png)

Fit a regression model.

##### French

``` r
m_acc_pred_error_gl_file <- file.path("..", "data", "model_fits", "m_acc_pred_error_Grandes_Lignes.rda")

if (file.exists(m_acc_pred_error_gl_file)) {
  load(m_acc_pred_error_gl_file)
} else {
  
  pred_gl_reg <- (
    pred_acc_error
    [course == "French"]
    [sample(.N, 1e6)]
    [, .(prediction_label, abs_acc_pred_error, user_id, fact_id)]
  )
  
  m_acc_pred_error_gl <- lmer(abs_acc_pred_error ~ prediction_label + 
                               (1 | user_id) + (1 | fact_id),
                             data = pred_gl_reg,
                             control = lmerControl(optimizer ="bobyqa"))
    
  save(m_acc_pred_error_gl, file = m_acc_pred_error_gl_file)
}

summary(m_acc_pred_error_gl)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs_acc_pred_error ~ prediction_label + (1 | user_id) + (1 |  
    ##     fact_id)
    ##    Data: pred_gl_reg
    ## Control: lmerControl(optimizer = "bobyqa")
    ## 
    ## REML criterion at convergence: -2363469
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -7.2812 -0.4159 -0.0378  0.4040  8.0239 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  user_id  (Intercept) 0.0004139 0.02034 
    ##  fact_id  (Intercept) 0.0004236 0.02058 
    ##  Residual             0.0051894 0.07204 
    ## Number of obs: 1000000, groups:  user_id, 40974; fact_id, 22910
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     4.900e-01  2.656e-04  6.217e+04 1844.99
    ## prediction_labelDomain         -1.095e-02  2.277e-04  9.749e+05  -48.09
    ## prediction_labelFact           -1.586e-02  2.292e-04  9.751e+05  -69.17
    ## prediction_labelLearner        -1.256e-02  2.310e-04  9.757e+05  -54.37
    ## prediction_labelFact & Learner -2.083e-02  2.321e-04  9.760e+05  -89.73
    ##                                Pr(>|t|)    
    ## (Intercept)                      <2e-16 ***
    ## prediction_labelDomain           <2e-16 ***
    ## prediction_labelFact             <2e-16 ***
    ## prediction_labelLearner          <2e-16 ***
    ## prediction_labelFact & Learner   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) prdc_D prdc_F prdc_L
    ## prdctn_lblD -0.429                     
    ## prdctn_lblF -0.420  0.497              
    ## prdctn_lblL -0.416  0.493  0.490       
    ## prdctn_lF&L -0.408  0.490  0.488  0.486

Compare different prediction types to each other:

``` r
ht_acc_gl <- glht(m_acc_pred_error_gl, linfct = mcp(prediction_label = "Tukey"))
summary(ht_acc_gl)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = abs_acc_pred_error ~ prediction_label + (1 | user_id) + 
    ##     (1 | fact_id), data = pred_gl_reg, control = lmerControl(optimizer = "bobyqa"))
    ## 
    ## Linear Hypotheses:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## Domain - Default == 0         -0.0109525  0.0002277 -48.092   <1e-10 ***
    ## Fact - Default == 0           -0.0158576  0.0002292 -69.173   <1e-10 ***
    ## Learner - Default == 0        -0.0125564  0.0002310 -54.367   <1e-10 ***
    ## Fact & Learner - Default == 0 -0.0208270  0.0002321 -89.730   <1e-10 ***
    ## Fact - Domain == 0            -0.0049051  0.0002293 -21.396   <1e-10 ***
    ## Learner - Domain == 0         -0.0016039  0.0002310  -6.944   <1e-10 ***
    ## Fact & Learner - Domain == 0  -0.0098745  0.0002322 -42.532   <1e-10 ***
    ## Learner - Fact == 0            0.0033012  0.0002324  14.202   <1e-10 ***
    ## Fact & Learner - Fact == 0    -0.0049694  0.0002334 -21.290   <1e-10 ***
    ## Fact & Learner - Learner == 0 -0.0082706  0.0002348 -35.222   <1e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Inspect the model's residuals:

``` r
qqnorm(resid(m_acc_pred_error_gl))
qqline(resid(m_acc_pred_error_gl), col = "red")
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-39-1.png)

``` r
plot(m_acc_pred_error_gl)
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-40-1.png)

##### English

``` r
m_acc_pred_error_ss_file <- file.path("..", "data", "model_fits", "m_acc_pred_error_Stepping_Stones.rda")

if (file.exists(m_acc_pred_error_ss_file)) {
  load(m_acc_pred_error_ss_file)
} else {
  
  pred_ss_reg <- (
    pred_acc_error
    [course == "English"]
    [sample(.N, 1e6)]
    [, .(prediction_label, abs_acc_pred_error, user_id, fact_id)]
  )
  
  m_acc_pred_error_ss <- lmer(abs_acc_pred_error ~ prediction_label + 
                               (1 | user_id) + (1 | fact_id),
                             data = pred_ss_reg,
                             control = lmerControl(optimizer ="bobyqa")
  )
  
  save(m_acc_pred_error_ss, file = m_acc_pred_error_ss_file)
}

summary(m_acc_pred_error_ss)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs_acc_pred_error ~ prediction_label + (1 | user_id) + (1 |  
    ##     fact_id)
    ##    Data: pred_ss_reg
    ## Control: lmerControl(optimizer = "bobyqa")
    ## 
    ## REML criterion at convergence: -2365380
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -7.2429 -0.3226 -0.0256  0.3277  9.3766 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  user_id  (Intercept) 0.0004519 0.02126 
    ##  fact_id  (Intercept) 0.0005107 0.02260 
    ##  Residual             0.0050246 0.07088 
    ## Number of obs: 1000000, groups:  user_id, 85899; fact_id, 45529
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     4.757e-01  2.298e-04  1.102e+05 2069.76
    ## prediction_labelDomain         -1.302e-02  2.281e-04  9.646e+05  -57.10
    ## prediction_labelFact           -1.546e-02  2.291e-04  9.640e+05  -67.46
    ## prediction_labelLearner        -1.432e-02  2.301e-04  9.635e+05  -62.21
    ## prediction_labelFact & Learner -1.969e-02  2.312e-04  9.630e+05  -85.17
    ##                                Pr(>|t|)    
    ## (Intercept)                      <2e-16 ***
    ## prediction_labelDomain           <2e-16 ***
    ## prediction_labelFact             <2e-16 ***
    ## prediction_labelLearner          <2e-16 ***
    ## prediction_labelFact & Learner   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) prdc_D prdc_F prdc_L
    ## prdctn_lblD -0.497                     
    ## prdctn_lblF -0.491  0.499              
    ## prdctn_lblL -0.489  0.497  0.494       
    ## prdctn_lF&L -0.483  0.494  0.492  0.490

Compare different prediction types to each other:

``` r
ht_acc_ss <- glht(m_acc_pred_error_ss, linfct = mcp(prediction_label = "Tukey"))
summary(ht_acc_ss)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = abs_acc_pred_error ~ prediction_label + (1 | user_id) + 
    ##     (1 | fact_id), data = pred_ss_reg, control = lmerControl(optimizer = "bobyqa"))
    ## 
    ## Linear Hypotheses:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## Domain - Default == 0         -0.0130225  0.0002281 -57.099   <1e-05 ***
    ## Fact - Default == 0           -0.0154554  0.0002291 -67.456   <1e-05 ***
    ## Learner - Default == 0        -0.0143165  0.0002301 -62.211   <1e-05 ***
    ## Fact & Learner - Default == 0 -0.0196905  0.0002312 -85.175   <1e-05 ***
    ## Fact - Domain == 0            -0.0024329  0.0002289 -10.630   <1e-05 ***
    ## Learner - Domain == 0         -0.0012941  0.0002299  -5.630   <1e-05 ***
    ## Fact & Learner - Domain == 0  -0.0066681  0.0002310 -28.868   <1e-05 ***
    ## Learner - Fact == 0            0.0011388  0.0002310   4.931   <1e-05 ***
    ## Fact & Learner - Fact == 0    -0.0042352  0.0002319 -18.264   <1e-05 ***
    ## Fact & Learner - Learner == 0 -0.0053740  0.0002329 -23.076   <1e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Inspect the model's residuals:

``` r
qqnorm(resid(m_acc_pred_error_ss))
qqline(resid(m_acc_pred_error_ss), col = "red")
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-43-1.png)

``` r
plot(m_acc_pred_error_ss)
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-44-1.png)

##### Comparison

``` r
ht_acc_gl_tidy <- broom::tidy(confint(ht_acc_gl))
ht_acc_ss_tidy <- broom::tidy(confint(ht_acc_ss))
setDT(ht_acc_gl_tidy)
setDT(ht_acc_ss_tidy)

ht_acc_both_tidy <- rbind(ht_acc_gl_tidy[, course := "French"],
                      ht_acc_ss_tidy[, course := "English"])
```

``` r
p_acc_pred_error_comp <- ggplot(ht_acc_both_tidy, aes(x = lhs, y = estimate, ymin = conf.low, ymax = conf.high, colour = course)) +
  geom_hline(yintercept = 0, linetype = "11", colour = "grey60") +
  geom_errorbar(width = 0.1) + 
  geom_point() +
  labs(x = "Linear hypotheses",
       y = "Estimate",
       caption = "Tukey's range test. Error bars show 95% family-wise confidence level.",
       colour = "Course") +
  coord_flip()

p_acc_pred_error_comp
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-46-1.png)

``` r
ggsave(plot = p_acc_pred_error_comp, file.path("..", "output", "acc_prediction_error_comparisons.png"),
       device = "png", width = 7.5, height = 5)

rm(p_acc_pred_error_comp)
```

##### Summary plot

``` r
pred_acc_error_avg[, prediction_rank := frank(-mae), by = .(course)]

annotation_df_ss <- data.table(
  course = rep("English", 10),
  start = c(1, 1, 1, 1,
            2, 2, 2,
            3, 3,
            4
  ),
  end = c(2, 3, 4, 5,
          3, 4, 5,
          4, 5,
          5
  ),
  y = seq(max(pred_acc_error_avg$mae)*1.01 + .01125, max(pred_acc_error_avg$mae)*1.01, by = -.00125),
  label = c("p < .001", "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001",
            "p < .001")
)

annotation_df_gl <- data.table(
  course = rep("French", 10),
  start = c(1, 1, 1, 1,
            2, 2, 2,
            3, 3,
            4
  ),
  end = c(2, 3, 4, 5,
          3, 4, 5,
          4, 5,
          5
  ),
  y = seq(max(pred_acc_error_avg$mae)*1.01 + .01125, max(pred_acc_error_avg$mae)*1.01, by = -.00125),
  label = c("p < .001", "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001",
            "p < .001")
)

annotation_df_acc <- rbind(annotation_df_ss, annotation_df_gl)
annotation_df_acc[, label := factor(label, levels = c("p < .001", "p < .01", "p < .05", "n.s."))]

p_acc_pred_error_summary <- ggplot(pred_acc_error_avg, aes(x = prediction_rank, y = mae)) +
  facet_grid(~ course) +
  geom_line(data = annotation_df_acc,
            aes(x = 1, y = .45, lty = label, alpha = label, colour = NULL)) + # Dummy line to get legend
  geom_line(aes(colour = course, group = course)) +
  geom_errorbar(aes(ymin = mae - ae_se, ymax = mae + ae_se), width = 0) +
  geom_point(aes(colour = course, group = course)) +
  geom_label(aes(label = prediction_label), 
             colour = "black", 
             alpha = .9,
             label.size = NA, 
             nudge_y = -.004) +
  labs(x = NULL,
       y = "Absolute prediction error:\nresponse accuracy",
       colour = "Course") +
  scale_x_continuous(expand = expansion(add = .75), breaks = NULL) +
  scale_colour_manual(values = dataset_colours) +
  scale_linetype_manual(values = c("p < .001" = 1,
                                   "p < .01" = 5,
                                   "p < .05" = 2,
                                   "n.s." = 3),
                        name = "Pairwise comparison:") +
  scale_alpha_manual(values = c("p < .001" = 1,
                                "p < .01" = .75,
                                "p < .05" = .5, 
                                "n.s." = .25), 
                     name = "Pairwise comparison:") +
  guides(colour = "none") +
  ggsignif::geom_signif(data = annotation_df_acc,
                        aes(xmin = start, xmax = end, annotations = "", 
                            y_position = y, lty = label, alpha = label),
                        tip_length = 0,
                        manual = TRUE)  +
  theme(legend.position = "bottom",
        legend.justification = "right")
```

    ## Warning: Ignoring unknown aesthetics: xmin, xmax, annotations, y_position

``` r
p_acc_pred_error_summary
```

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-47-1.png)

``` r
ggsave(file.path("..", "output", "acc_absolute_prediction_error_summary.png"),
       device = "png", width = 10, height = 4)
```

##### Improvement

How big was the improvement from worst to best prediction method?

French:

``` r
# Absolute change
ht_acc_gl_tidy[lhs == "Fact & Learner - Default", estimate[1]]
```

    ## [1] -0.02082699

``` r
# % change
scales::percent(
  ht_acc_gl_tidy[lhs == "Fact & Learner - Default", estimate[1]] / fixef(m_acc_pred_error_gl)[[1]],
  accuracy = .1)
```

    ## [1] "-4.3%"

English:

``` r
# Absolute change
ht_acc_ss_tidy[lhs == "Fact & Learner - Default", estimate[1]]
```

    ## [1] -0.01969055

``` r
# % change
scales::percent(
  ht_acc_ss_tidy[lhs == "Fact & Learner - Default", estimate[1]] / fixef(m_acc_pred_error_ss)[[1]],
  accuracy = .1)
```

    ## [1] "-4.1%"

# Combined plot

``` r
(p_acc_pred_error_summary + p_rt_pred_error_summary) + 
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
```

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

![](05_evaluate_behavioural_predictions_files/figure-markdown_github/unnamed-chunk-50-1.png)

``` r
ggsave(file.path("..", "output", "beh_absolute_prediction_error_summary.png"),
       device = "png", width = 10, height = 8)
```

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

    ## Warning in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, :
    ## 'big.mark' and 'decimal.mark' are both '.', which could be confusing

# Session info

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] dplyr_1.0.7       multcomp_1.4-10   TH.data_1.0-10   
    ##  [4] MASS_7.3-51.4     survival_2.44-1.1 mvtnorm_1.1-1    
    ##  [7] lmerTest_3.1-0    lme4_1.1-21       Matrix_1.2-18    
    ## [10] wesanderson_0.3.6 patchwork_1.1.1   ggplot2_3.3.5    
    ## [13] stringr_1.4.0     furrr_0.1.0       future_1.13.0    
    ## [16] purrr_0.3.2       tidyr_1.0.0       data.table_1.13.6
    ## [19] fst_0.9.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] zoo_1.8-6           tidyselect_1.1.1    xfun_0.21          
    ##  [4] listenv_0.7.0       splines_3.6.3       lattice_0.20-41    
    ##  [7] colorspace_1.4-1    vctrs_0.3.8         generics_0.1.0     
    ## [10] htmltools_0.3.6     mgcv_1.8-28         yaml_2.2.0         
    ## [13] utf8_1.1.4          rlang_0.4.10        pillar_1.6.3       
    ## [16] nloptr_1.2.1        glue_1.4.2          withr_2.3.0        
    ## [19] DBI_1.1.0           lifecycle_1.0.1     ggsignif_0.5.0     
    ## [22] munsell_0.5.0       gtable_0.3.0        codetools_0.2-16   
    ## [25] evaluate_0.14       labeling_0.3        knitr_1.23         
    ## [28] parallel_3.6.3      fansi_0.4.0         broom_0.5.2        
    ## [31] Rcpp_1.0.6          backports_1.1.4     scales_1.1.1       
    ## [34] jsonlite_1.6        farver_2.1.0        digest_0.6.19      
    ## [37] stringi_1.4.3       numDeriv_2016.8-1.1 grid_3.6.3         
    ## [40] tools_3.6.3         sandwich_2.5-1      magrittr_2.0.1     
    ## [43] tibble_2.1.3        crayon_1.4.1        pkgconfig_2.0.2    
    ## [46] ellipsis_0.3.2      minqa_1.2.4         rmarkdown_2.6      
    ## [49] R6_2.4.0            globals_0.12.4      boot_1.3-25        
    ## [52] nlme_3.1-149        compiler_3.6.3
