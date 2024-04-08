Evaluate rate of forgetting predictions
================
Maarten van der Velde
Last updated: 2024-01-10

- [Overview](#overview)
- [Setup](#setup)
  - [Helper functions](#helper-functions)
- [Rate of forgetting](#rate-of-forgetting)
  - [Predicted rate of forgetting](#predicted-rate-of-forgetting)
    - [Distribution of predictions](#distribution-of-predictions)
    - [Predicted vs. observed values](#predicted-vs-observed-values)
    - [Prediction error](#prediction-error)
    - [Visualise prediction error](#visualise-prediction-error)
- [Session info](#session-info)

# Overview

This notebook evaluates the predicted rates of forgetting. The
predictions are compared to the value derived at the end of each
learning session to determine their accuracy.

# Setup

``` r
library(fst)
library(data.table)
library(tidyr)
library(purrr)
library(furrr)
library(stringr)
library(ggplot2)
library(wesanderson)
library(lme4)
library(lmerTest)
library(multcomp)
```

``` r
source(file.path("..", "scripts", "99_slimstampen_model_funs.R"))
```

``` r
future::plan("multisession", workers = 6) # Set to desired number of cores
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
load_predictions <- function(course) {
  
  pred_user <- read_fst(file.path("..", "data", "predictions", paste0("pred_v_obs_user_", str_replace_all(course, " ", "_"), ".fst")))
  pred_fact <- read_fst(file.path("..", "data", "predictions", paste0("pred_v_obs_fact_", str_replace_all(course, " ", "_"), ".fst")))
  pred_fact_user <- read_fst(file.path("..", "data", "predictions", paste0("pred_fact_and_user_", str_replace_all(course, " ", "_"), ".fst")))
  setDT(pred_user)
  setDT(pred_fact)
  setDT(pred_fact_user)
  
  pred_domain <- mean(unique(pred_fact, by = c("fact_id"))$pred_fact)
  pred_default <- 0.3
  
  # Combine
  pred_all <- merge(pred_user, pred_fact, by = c("user_id", "fact_id", "alpha", "n_reps"), all = TRUE)
  pred_all <- merge(pred_all, pred_fact_user, by = c("user_id", "fact_id", "alpha"), all = TRUE)
  pred_all[, pred_default := pred_default]
  pred_all[, pred_domain := pred_domain]
  
  pred_obs_long <- pivot_longer(pred_all, 
                                cols = pred_user:pred_domain,
                                names_to = "prediction_type",
                                names_prefix = "pred\\_")
  
  setDT(pred_obs_long)
  
  # Remove NA predictions and predictions without corresponding observations
  pred_obs_long <- pred_obs_long[!is.na(value)]
  pred_obs_long <- pred_obs_long[!is.na(alpha)]
  
  # Remove duplicates
  pred_obs_long <- unique(pred_obs_long)
  
  # Set proper labels
  condition_labels <- data.table(prediction_type = c("default", "domain", "fact", "user", "fact_user"),
                                 prediction_label = factor(c("Default", "Domain", "Fact", "Learner", "Fact & Learner"),
                                                           levels = c("Default", "Domain", "Fact", "Learner", "Fact & Learner")))
  pred_obs_long <- pred_obs_long[condition_labels, on = .(prediction_type)]
  
  return(pred_obs_long)
}
```

# Rate of forgetting

## Predicted rate of forgetting

``` r
pred_gl <- load_predictions("Grandes Lignes")
pred_ss <- load_predictions("Stepping Stones")

pred_both <- rbind(pred_gl[, course := "French"],
                   pred_ss[, course := "English"])
```

### Distribution of predictions

``` r
p_rof_dist <- ggplot(pred_both, aes(x = value, fill = prediction_label)) +
  facet_grid(prediction_label ~ course, scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .01) +
  guides(fill = "none") +
  labs(x = "Predicted rate-of-forgetting",
       y = "Density") +
  scale_fill_manual(values = condition_colours)

p_rof_dist
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave(plot = p_rof_dist, file.path("..", "output", "rof_predictions_distribution.png"),
       device = "png", width = 5, height = 7.5)
```

### Predicted vs. observed values

We compare the predicted rate of forgetting to the “observed” rate of
forgetting, i.e., the rate of forgetting that was estimated at the end
of the learning sequence.

To assess the accuracy of predictions, we compute the mean absolute
error (MAE) as an aggregate statistic, as well as the absolute error
(AE) of each individual prediction.

``` r
pred_mae <- pred_both[, .(mae = mean(abs(alpha - value)),
                          ae_se = sd(abs(alpha - value))/.N), 
                      by = .(course, prediction_label)]

  
n_obs <- pred_both[, .N, by = .(course, prediction_label)]
```

Plot predicted vs. observed values:

``` r
rof_min <- 0
rof_max <- 1
rof_breaks <- seq(0.1, 0.9, by = .2)

p_rof_pred_v_obs <- ggplot(pred_both,
                         aes(x = value, y = alpha, colour = prediction_label)) +
    facet_grid(course ~ prediction_label) +
    geom_hline(yintercept = 0.3, lty = 2) +
    geom_vline(xintercept = 0.3, lty = 2) +
    geom_abline(slope = 1, intercept = 0, lty = 3, alpha = 0.75) +
    geom_point(alpha = .1, size = .1, pch = ".") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
  geom_label(data = pred_mae,
             aes(label = paste("MAE =", formatC(mae, digits = 3, flag = "#"))),
             x = rof_max, y = rof_min,
             hjust = 1, colour = "NA", size = 3,
             alpha = .9,
             label.size = NA) +
  geom_text(data = pred_mae,
            aes(label = paste("MAE =", formatC(mae, digits = 3, flag = "#"))),
            x = rof_max, y = rof_min,
            hjust = 1, colour = "black", size = 3) +
  geom_label(data = n_obs,
             aes(label = paste("n =", scales::comma(N))),
             x = rof_max,
             y = rof_max,
             hjust = 1, colour = "NA", size = 3,
             alpha = .9,
             label.size = NA) +
  geom_text(data = n_obs,
            aes(label = paste("n =", scales::comma(N))),
            x = rof_max,
            y = rof_max,
            hjust = 1, colour = "black", size = 3) +
  guides(colour = "none") +
  labs(x = "Predicted rate-of-forgetting α",
       y = "Observed rate-of-forgetting α") +
  coord_fixed(ratio = 1, xlim = c(rof_min, rof_max), ylim = c(rof_min, rof_max)) +
  scale_x_continuous(breaks = rof_breaks) +
  scale_y_continuous(breaks = rof_breaks) +
  scale_colour_manual(values = condition_colours)

p_rof_pred_v_obs
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave(plot = p_rof_pred_v_obs, file.path("..", "output", "rof_predicted_vs_observed.png"),
  device = "png", width = 10, height = 4.5)

rm(p_rof_pred_v_obs)
```

### Prediction error

Calculate prediction error:

``` r
pred_both[, prediction_error := value - alpha]
```

Distribution of prediction error:

``` r
p_rof_pred_error <- ggplot(pred_both, aes(x = prediction_error, fill = prediction_label)) +
  facet_grid(prediction_label ~ course, scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .01) +
  guides(fill = "none") +
  labs(x = "Rate-of-forgetting prediction error (predicted - observed)",
       y = "Density") +
  scale_fill_manual(values = condition_colours)

p_rof_pred_error
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave(plot = p_rof_pred_error, file.path("..", "output", "rof_prediction_error.png"),
       device = "png", width = 5, height = 7.5)
```

#### Absolute prediction error

To compare the magnitude of prediction errors between prediction
methods, we look at absolute prediction error.

``` r
pred_both[, abs_prediction_error := abs(prediction_error)]
```

Distribution of absolute prediction error:

``` r
p_rof_abs_pred_error <- ggplot(pred_both, aes(x = abs_prediction_error, fill = prediction_label)) +
  facet_grid(prediction_label ~ course, scales = "free_y") +
  geom_histogram(aes(y = ..density..), binwidth = .01) +
  guides(fill = "none") +
  labs(x = "Absolute rate-of-forgetting prediction error",
       y = "Density") +
  scale_fill_manual(values = condition_colours)

p_rof_abs_pred_error
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggsave(plot = p_rof_abs_pred_error, file.path("..", "output", "rof_absolute_prediction_error.png"),
       device = "png", width = 5, height = 7.5)
```

``` r
pred_error_summarised <- pred_both[, .(error_mean = mean(abs_prediction_error), error_se = sd(abs_prediction_error)/.N), by = .(course, prediction_label)]

ggplot(pred_error_summarised, aes(x = prediction_label, y = error_mean, colour = course)) +
  geom_boxplot(data = pred_both,
               aes(y = abs_prediction_error, group = interaction(course, prediction_label)),
               colour = "grey70",
               width = .25,
               outlier.shape = NA,
               position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = error_mean - error_se, ymax = error_mean + error_se), width = 0, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  coord_cartesian(ylim = c(0, 0.175)) +
  labs(x = "Method",
       y = "Absolute rate-of-forgetting prediction error",
       colour = "Course")
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Fit a regression model on absolute prediction error. The whole data set
is too big to fit in a reasonable time, so we fit the model to a random
subset of 1M predictions (which already takes ~24hrs).

##### Grandes Lignes

``` r
m_rof_pred_error_gl_file <- file.path("..", "data", "model_fits", "m_pred_error_Grandes_Lignes_1e6.rda")

if (file.exists(m_rof_pred_error_gl_file)) {
  load(m_rof_pred_error_gl_file)
} else {
  
  pred_gl_reg <- (
    pred_both
    [course == "French"]
    [sample(.N, 1e6), .(prediction_label, abs_prediction_error, user_id, fact_id)]
  )
  
  m_rof_pred_error_gl <- lmer(abs_prediction_error ~ prediction_label + 
                                (1 | user_id) + (1 | fact_id),
                              data = pred_gl_reg)
  
  save(m_rof_pred_error_gl, file = m_rof_pred_error_gl_file)
}

summary(m_rof_pred_error_gl)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs_prediction_error ~ prediction_label + (1 | user_id) + (1 |  
    ##     fact_id)
    ##    Data: pred_gl_reg
    ## 
    ## REML criterion at convergence: -3273096
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6381 -0.5902 -0.1803  0.3945 14.5510 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  user_id  (Intercept) 0.0001461 0.01209 
    ##  fact_id  (Intercept) 0.0001916 0.01384 
    ##  Residual             0.0020919 0.04574 
    ## Number of obs: 1000000, groups:  user_id, 40965; fact_id, 22884
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     6.252e-02  1.701e-04  5.391e+04 367.666
    ## prediction_labelDomain         -1.430e-03  1.446e-04  9.734e+05  -9.892
    ## prediction_labelFact           -1.256e-02  1.454e-04  9.736e+05 -86.340
    ## prediction_labelLearner        -4.216e-03  1.465e-04  9.744e+05 -28.786
    ## prediction_labelFact & Learner -1.219e-02  1.474e-04  9.746e+05 -82.745
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
    ## prdctn_lblD -0.425                     
    ## prdctn_lblF -0.416  0.496              
    ## prdctn_lblL -0.413  0.493  0.490       
    ## prdctn_lF&L -0.404  0.490  0.488  0.486

Compare different prediction types to each other:

``` r
ht_gl <- glht(m_rof_pred_error_gl, linfct = mcp(prediction_label = "Tukey"))
summary(ht_gl)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = abs_prediction_error ~ prediction_label + (1 | 
    ##     user_id) + (1 | fact_id), data = pred_gl_reg)
    ## 
    ## Linear Hypotheses:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## Domain - Default == 0         -0.0014305  0.0001446  -9.892   <0.001 ***
    ## Fact - Default == 0           -0.0125572  0.0001454 -86.340   <0.001 ***
    ## Learner - Default == 0        -0.0042160  0.0001465 -28.786   <0.001 ***
    ## Fact & Learner - Default == 0 -0.0121938  0.0001474 -82.745   <0.001 ***
    ## Fact - Domain == 0            -0.0111266  0.0001456 -76.442   <0.001 ***
    ## Learner - Domain == 0         -0.0027855  0.0001466 -18.999   <0.001 ***
    ## Fact & Learner - Domain == 0  -0.0107633  0.0001475 -72.977   <0.001 ***
    ## Learner - Fact == 0            0.0083412  0.0001474  56.590   <0.001 ***
    ## Fact & Learner - Fact == 0     0.0003633  0.0001482   2.452    0.102    
    ## Fact & Learner - Learner == 0 -0.0079779  0.0001490 -53.534   <0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Inspect the model’s residuals:

``` r
qqnorm(resid(m_rof_pred_error_gl))
qqline(resid(m_rof_pred_error_gl), col = "red")
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plot(m_rof_pred_error_gl)
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

The QQ plot indicates quite a strong skew, which is not surprising,
given that the distribution of absolute error is bounded by zero on the
left but unbounded on the right. Assuming a Gamma distribution may be
better, but models that use a Gamma distribution do not converge here.
The LMER also gives a sufficiently accurate estimate of the means.

##### Stepping Stones

``` r
m_rof_pred_error_ss_file <- file.path("..", "data", "model_fits", "m_pred_error_Stepping_Stones_1e6.rda")

if (file.exists(m_rof_pred_error_ss_file)) {
  load(m_rof_pred_error_ss_file)
} else {
  
  pred_ss_reg <- (
    pred_both
    [course == "English"]
    [sample(.N, 1e6), .(prediction_label, abs_prediction_error, user_id, fact_id)]
  )
  
  m_pred_error <- lmer(abs_prediction_error ~ prediction_label + 
                         (1 | user_id) + (1 | fact_id),
                       data = pred_ss_reg)
  
  save(m_pred_error, file = m_rof_pred_error_ss_file)
}

summary(m_pred_error)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: abs_prediction_error ~ prediction_label + (1 | user_id) + (1 |  
    ##     fact_id)
    ##    Data: pred_ss_reg
    ## 
    ## REML criterion at convergence: -3582447
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1157 -0.5163 -0.1535  0.2167 16.6397 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  user_id  (Intercept) 0.0001343 0.01159 
    ##  fact_id  (Intercept) 0.0001022 0.01011 
    ##  Residual             0.0014998 0.03873 
    ## Number of obs: 1000000, groups:  user_id, 86084; fact_id, 45580
    ## 
    ## Fixed effects:
    ##                                  Estimate Std. Error         df t value
    ## (Intercept)                     5.297e-02  1.188e-04  1.484e+05  445.85
    ## prediction_labelDomain         -2.163e-03  1.244e-04  9.714e+05  -17.39
    ## prediction_labelFact           -7.743e-03  1.249e-04  9.710e+05  -62.02
    ## prediction_labelLearner        -5.282e-03  1.255e-04  9.705e+05  -42.08
    ## prediction_labelFact & Learner -7.926e-03  1.260e-04  9.701e+05  -62.91
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
    ## prdctn_lblD -0.522                     
    ## prdctn_lblF -0.517  0.497              
    ## prdctn_lblL -0.513  0.495  0.493       
    ## prdctn_lF&L -0.509  0.493  0.492  0.489

Compare different prediction types to each other:

``` r
ht_ss <- glht(m_pred_error, linfct = mcp(prediction_label = "Tukey"))
summary(ht_ss)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = abs_prediction_error ~ prediction_label + (1 | 
    ##     user_id) + (1 | fact_id), data = pred_ss_reg)
    ## 
    ## Linear Hypotheses:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## Domain - Default == 0         -0.0021631  0.0001244 -17.393   <1e-04 ***
    ## Fact - Default == 0           -0.0077430  0.0001249 -62.018   <1e-04 ***
    ## Learner - Default == 0        -0.0052823  0.0001255 -42.083   <1e-04 ***
    ## Fact & Learner - Default == 0 -0.0079259  0.0001260 -62.907   <1e-04 ***
    ## Fact - Domain == 0            -0.0055799  0.0001249 -44.665   <1e-04 ***
    ## Learner - Domain == 0         -0.0031192  0.0001256 -24.834   <1e-04 ***
    ## Fact & Learner - Domain == 0  -0.0057629  0.0001261 -45.717   <1e-04 ***
    ## Learner - Fact == 0            0.0024607  0.0001261  19.518   <1e-04 ***
    ## Fact & Learner - Fact == 0    -0.0001830  0.0001265  -1.447    0.597    
    ## Fact & Learner - Learner == 0 -0.0026437  0.0001271 -20.798   <1e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

Inspect the model’s residuals:

``` r
qqnorm(resid(m_pred_error))
qqline(resid(m_pred_error), col = "red")
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

##### Comparison

``` r
ht_gl_tidy <- broom::tidy(confint(ht_gl))
ht_ss_tidy <- broom::tidy(confint(ht_ss))
setDT(ht_gl_tidy)
setDT(ht_ss_tidy)

ht_both_tidy <- rbind(ht_gl_tidy[, course := "French"],
                      ht_ss_tidy[, course := "English"])
```

``` r
p_rof_pred_error_comp <- ggplot(ht_both_tidy, aes(x = contrast, y = estimate, ymin = conf.low, ymax = conf.high, colour = course)) +
  geom_hline(yintercept = 0, linetype = "11", colour = "grey60") +
  geom_errorbar(width = 0.1) + 
  geom_point() +
  labs(x = "Linear hypotheses",
       y = "Estimate",
       caption = "Tukey's range test. Error bars show 95% family-wise confidence level.",
       colour = "Course") +
    coord_flip()

p_rof_pred_error_comp
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
ggsave(plot = p_rof_pred_error_comp, file.path("..", "output", "rof_prediction_error_comparisons.png"),
       device = "png", width = 7.5, height = 5)
```

##### Summary plot

``` r
# Set significance level of comparisons manually, based on model output
pred_error_summarised$comparison <- c("***")
pred_error_summarised[c(3, 8), comparison := NA]
pred_error_summarised[c(5, 10), comparison := "n.s."]

# Add fitted values
pred_error_summarised[course == "French", error_fitted := predict(m_rof_pred_error_gl,
                                                newdata = pred_error_summarised[course == "French"],
                                                re.form = NA, 
                                                type = "response")]
pred_error_summarised[course == "English", error_fitted := predict(m_pred_error,
                                                newdata = pred_error_summarised[course == "English"],
                                                re.form = NA, 
                                                type = "response")]

p_rof_abs_pred_error_summ <- ggplot(pred_error_summarised, aes(x = reorder(prediction_label, -error_mean), y = error_mean, colour = course)) +
  geom_errorbar(aes(ymin = error_mean - error_se, ymax = error_mean + error_se), width = 0) +
  geom_line(aes(group = course), lty = 2) +
  geom_point() +
  geom_text(aes(label = comparison), 
            colour = "black",
            position = position_nudge(x = .5, y = c(rep(0, 9), .001)),
            hjust = .5) +
  labs(x = "Method",
       y = "Absolute rate-of-forgetting prediction error",
       colour = "Course") +
  scale_colour_manual(values = dataset_colours) +
  theme(legend.position = c(.85, .85))

p_rof_abs_pred_error_summ
```

    ## Warning: Removed 2 rows containing missing values (`geom_text()`).

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
ggsave(plot = p_rof_abs_pred_error_summ, file.path("..", "output", "rof_absolute_prediction_error_summary.png"),
       device = "png", width = 7.5, height = 4.5)
```

    ## Warning: Removed 2 rows containing missing values (`geom_text()`).

``` r
pred_error_summarised[, prediction_rank := frank(-error_mean), by = .(course)]

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
  y = seq(max(pred_error_summarised$error_mean)*1.01 + .00675, max(pred_error_summarised$error_mean)*1.01, by = -.00075),
  label = c("p < .001", "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001",
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
  y = seq(max(pred_error_summarised$error_mean)*1.01 + .00675, max(pred_error_summarised$error_mean)*1.01, by = -.00075),
  label = c("p < .001", "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001", "p < .001",
            "p < .001", "p < .001",
            "n.s.")
)

annotation_df_rof <- rbind(annotation_df_ss, annotation_df_gl)
annotation_df_rof[, label := factor(label, levels = c("p < .001", "p < .01", "p < .05", "n.s."))]

p_rof_pred_error_summary <- ggplot(pred_error_summarised, aes(x = prediction_rank, y = error_mean)) +
  facet_grid(~ course) +
  geom_line(data = annotation_df_rof,
            aes(x = 1, y = .05, lty = label, alpha = label, colour = NULL)) + # Dummy line to get legend
  geom_line(aes(colour = course, group = course)) +
  geom_errorbar(aes(ymin = error_mean - error_se, ymax = error_mean + error_se), width = 0) +
  geom_point(aes(colour = course, group = course)) +
  geom_label(aes(label = prediction_label), 
             colour = "black", 
             alpha = .9,
             label.size = NA, 
             nudge_y = -.0025) +
  labs(x = NULL,
       y = "Absolute prediction error:\nrate-of-forgetting α",
       colour = "Course") +
  scale_x_continuous(expand = expansion(add = .75), breaks = NULL) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = dataset_colours) +
  scale_linetype_manual(values = c("p < .001" = 1,
                                   "p < .01" = 5,
                                   "p < .05" = 2,
                                   "n.s." = 3),
                        drop = FALSE,
                        name = "Pairwise comparison:") +
  scale_alpha_manual(values = c("p < .001" = 1,
                                "p < .01" = .75,
                                "p < .05" = .5, 
                                "n.s." = .25),
                     drop = FALSE,
                     name = "Pairwise comparison:") +
  guides(colour = "none") +
  ggsignif::geom_signif(data = annotation_df_rof,
                        aes(xmin = start, xmax = end, annotations = "", 
                            y_position = y, lty = label, alpha = label),
                        tip_length = 0,
                        manual = TRUE)  +
  theme(legend.position = "bottom",
        legend.justification = "right")
```

    ## Warning in ggsignif::geom_signif(data = annotation_df_rof, aes(xmin = start, :
    ## Ignoring unknown aesthetics: xmin, xmax, annotations, and y_position

``` r
p_rof_pred_error_summary
```

    ## Warning: The following aesthetics were dropped during statistical transformation: xmin, xmax, y_position
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor?

    ## Warning: The following aesthetics were dropped during statistical transformation: xmin, xmax, y_position
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor?

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_absolute_prediction_error_summary.png"),
       device = "png", width = 10, height = 8)
```

    ## Warning: The following aesthetics were dropped during statistical transformation: xmin, xmax, y_position
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor?
    ## The following aesthetics were dropped during statistical transformation: xmin, xmax, y_position
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor?

#### Improvement

How big was the improvement from worst to best prediction method?

French:

``` r
# Absolute change
ht_gl_tidy[contrast == "Fact - Default", estimate[1]]
```

    ## [1] -0.01255716

``` r
# % change
scales::percent(
  ht_gl_tidy[contrast == "Fact - Default", estimate[1]] / fixef(m_rof_pred_error_gl)[[1]],
  accuracy = .1)
```

    ## [1] "-20.1%"

English:

``` r
# Absolute change
ht_ss_tidy[contrast == "Fact & Learner - Default", estimate[1]]
```

    ## [1] -0.007925937

``` r
# % change
scales::percent(
  ht_ss_tidy[contrast == "Fact & Learner - Default", estimate[1]] / fixef(m_pred_error)[[1]],
  accuracy = .1)
```

    ## [1] "-15.0%"

### Visualise prediction error

#### By learner

``` r
user_freq <- pred_both[, .N, by = .(course, prediction_label, user_id)]

pred_user_freq <- pred_both[user_freq[N > 50], on = .(course, prediction_label, user_id)]

pred_user_q <- pred_user_freq[, .(stat = c("whisker_low", "q25", "median", "q75", "whisker_high"),
                                  value = boxplot.stats(prediction_error, do.conf = FALSE, do.out = FALSE)$stats), by = .(course, prediction_label, user_id)]

pred_user_q <- pivot_wider(pred_user_q, names_from = "stat", values_from = "value")

pred_user_q <- pred_user_q %>%
  arrange(course, prediction_label, median) %>%
  group_by(course, prediction_label) %>%
  mutate(user_order = (1:n())/n())
```

``` r
ggplot(pred_user_q, aes(x = user_order)) +
  facet_grid(course ~ prediction_label) +
  geom_ribbon(aes(ymin = whisker_low, ymax = whisker_high, fill = course), alpha = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = course), alpha = .5) +
  geom_line(aes(y = median), lwd = 1) +
  geom_hline(data = NULL, yintercept = 0, lty = 3) +
    labs(x = "Learners",
       y = "Rate-of-forgetting prediction error\n(predicted - observed)") +
  scale_fill_manual(values = dataset_colours) +
  guides(fill = "none") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_prediction_error_by_learner.png"),
  device = "png", width = 10, height = 4.5)
```

Absolute error:

``` r
pred_user_q <- pred_user_freq[, .(stat = c("whisker_low", "q25", "median", "q75", "whisker_high"),
                                  value = boxplot.stats(abs_prediction_error, do.conf = FALSE, do.out = FALSE)$stats), by = .(course, prediction_label, user_id)]

pred_user_q <- pivot_wider(pred_user_q, names_from = "stat", values_from = "value")

pred_user_q <- pred_user_q %>%
  arrange(course, prediction_label, median) %>%
  group_by(course, prediction_label) %>%
  mutate(user_order = (1:n())/n())

ggplot(pred_user_q, aes(x = user_order)) +
  facet_grid(course ~ prediction_label) +
  geom_ribbon(aes(ymin = whisker_low, ymax = whisker_high, fill = course), alpha = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = course), alpha = .5) +
  geom_line(aes(y = median), lwd = 1) +
    labs(x = "Learners",
       y = "Absolute rate-of-forgetting prediction error") +
  scale_fill_manual(values = dataset_colours) +
  guides(fill = "none") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_abs_prediction_error_by_learner.png"),
  device = "png", width = 10, height = 4.5)
```

``` r
ggplot(pred_user_q, aes(x = user_order, group = prediction_label, colour = prediction_label)) +
  facet_grid(~ course) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = prediction_label), colour = NA, alpha = .15) +
  geom_line(aes(y = median), lwd = 1) +
    labs(x = "Learners",
       y = "Absolute rate-of-forgetting prediction error",
       colour = "Prediction\nmethod",
       fill = "Prediction\nmethod") +
  scale_colour_manual(values = condition_colours) +
  scale_fill_manual(values = condition_colours) +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_abs_prediction_error_by_learner_condensed.png"),
  device = "png", width = 10, height = 4.5)
```

#### By fact

``` r
fact_freq <- pred_both[, .N, by = .(course, prediction_label, fact_id)]

pred_fact_freq <- pred_both[fact_freq[N > 50], on = .(course, prediction_label, fact_id)]

pred_fact_q <- pred_fact_freq[, .(stat = c("whisker_low", "q25", "median", "q75", "whisker_high"),
                                  value = boxplot.stats(prediction_error, do.conf = FALSE, do.out = FALSE)$stats), by = .(course, prediction_label, fact_id)]

pred_fact_q <- pivot_wider(pred_fact_q, names_from = "stat", values_from = "value")

pred_fact_q <- pred_fact_q %>%
  arrange(course, prediction_label, median) %>%
  group_by(course, prediction_label) %>%
  mutate(fact_order = (1:n())/n())
```

``` r
ggplot(pred_fact_q, aes(x = fact_order)) +
  facet_grid(course ~ prediction_label) +
  geom_ribbon(aes(ymin = whisker_low, ymax = whisker_high, fill = course), alpha = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = course), alpha = .5) +
  geom_line(aes(y = median), lwd = 1) +
  geom_hline(data = NULL, yintercept = 0, lty = 3) +
    labs(x = "Facts",
       y = "Rate-of-forgetting prediction error\n(predicted - observed)") +
  scale_fill_manual(values = dataset_colours) +
  guides(fill = "none") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_prediction_error_by_fact.png"),
  device = "png", width = 10, height = 4.5)
```

Absolute error:

``` r
pred_fact_q <- pred_fact_freq[, .(stat = c("whisker_low", "q25", "median", "q75", "whisker_high"),
                                  value = boxplot.stats(abs_prediction_error, do.conf = FALSE, do.out = FALSE)$stats), by = .(course, prediction_label, fact_id)]

pred_fact_q <- pivot_wider(pred_fact_q, names_from = "stat", values_from = "value")

pred_fact_q <- pred_fact_q %>%
  arrange(course, prediction_label, median) %>%
  group_by(course, prediction_label) %>%
  mutate(fact_order = (1:n())/n())

ggplot(pred_fact_q, aes(x = fact_order)) +
  facet_grid(course ~ prediction_label) +
  geom_ribbon(aes(ymin = whisker_low, ymax = whisker_high, fill = course), alpha = .3) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = course), alpha = .5) +
  geom_line(aes(y = median), lwd = 1) +
    labs(x = "Facts",
       y = "Absolute rate-of-forgetting prediction error") +
  scale_fill_manual(values = dataset_colours) +
  guides(fill = "none") +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_abs_prediction_error_by_fact.png"),
  device = "png", width = 10, height = 4.5)
```

``` r
ggplot(pred_fact_q, aes(x = fact_order, group = prediction_label, colour = prediction_label)) +
  facet_grid(~ course) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = prediction_label), colour = NA, alpha = .15) +
  geom_line(aes(y = median), lwd = 1) +
    labs(x = "Facts",
       y = "Absolute rate-of-forgetting prediction error",
       colour = "Prediction\nmethod",
       fill = "Prediction\nmethod") +
  scale_colour_manual(values = condition_colours) +
  scale_fill_manual(values = condition_colours) +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_abs_prediction_error_by_fact_condensed.png"),
  device = "png", width = 10, height = 4.5)
```

``` r
ggplot(pred_both, aes(x = abs_prediction_error, colour = prediction_label, fill = prediction_label)) +
  facet_grid(~ course) +
  geom_density(alpha = .1) +
      labs(x = "Absolute rate-of-forgetting prediction error",
       y = "Density",
       colour = "Prediction\nmethod",
       fill = "Prediction\nmethod") +
  scale_colour_manual(values = condition_colours) +
  scale_fill_manual(values = condition_colours) +
  coord_cartesian(ylim = c(0, 100), xlim = c(0, .25)) +
  NULL
```

![](04_evaluate_ROF_predictions_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
ggsave(file.path("..", "output", "rof_abs_prediction_error_density.png"),
  device = "png", width = 10, height = 4.5)
```

# Session info

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Sonoma 14.2.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/Amsterdam
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] fstcore_0.9.14    dplyr_1.1.3       multcomp_1.4-25   TH.data_1.1-2    
    ##  [5] MASS_7.3-60       survival_3.5-7    mvtnorm_1.2-3     lmerTest_3.1-3   
    ##  [9] lme4_1.1-34       Matrix_1.6-1.1    wesanderson_0.3.6 ggplot2_3.4.3    
    ## [13] stringr_1.5.0     furrr_0.3.1       future_1.33.0     purrr_1.0.2      
    ## [17] tidyr_1.3.0       data.table_1.14.8 fst_0.9.8        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.4        xfun_0.40           bslib_0.5.1        
    ##  [4] lattice_0.21-9      numDeriv_2016.8-1.1 vctrs_0.6.3        
    ##  [7] tools_4.3.1         generics_0.1.3      parallel_4.3.1     
    ## [10] sandwich_3.0-2      tibble_3.2.1        fansi_1.0.4        
    ## [13] pkgconfig_2.0.3     lifecycle_1.0.3     compiler_4.3.1     
    ## [16] farver_2.1.1        munsell_0.5.0       codetools_0.2-19   
    ## [19] htmltools_0.5.6     sass_0.4.7          yaml_2.3.7         
    ## [22] crayon_1.5.2        pillar_1.9.0        nloptr_2.0.3       
    ## [25] jquerylib_0.1.4     cachem_1.0.8        boot_1.3-28.1      
    ## [28] nlme_3.1-163        parallelly_1.36.0   tidyselect_1.2.0   
    ## [31] digest_0.6.33       stringi_1.7.12      listenv_0.9.0      
    ## [34] labeling_0.4.3      splines_4.3.1       fastmap_1.1.1      
    ## [37] grid_4.3.1          colorspace_2.1-0    cli_3.6.1          
    ## [40] magrittr_2.0.3      utf8_1.2.3          broom_1.0.5        
    ## [43] withr_2.5.1         backports_1.4.1     scales_1.2.1       
    ## [46] rmarkdown_2.25      globals_0.16.2      ggsignif_0.6.4     
    ## [49] zoo_1.8-12          evaluate_0.22       knitr_1.44         
    ## [52] mgcv_1.9-0          rlang_1.1.1         Rcpp_1.0.11        
    ## [55] glue_1.6.2          rstudioapi_0.15.0   minqa_1.2.6        
    ## [58] jsonlite_1.8.7      R6_2.5.1
