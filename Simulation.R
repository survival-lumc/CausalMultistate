### Loading libraries and functions --------------------------------------------
### ----------------------------------------------------------------------------
library(survival)
library(mstate)
library(tidyverse)
library(matrixStats) 

source("SimulationFunctions.R")


### Set up parameters ----------------------------------------------------------
### ----------------------------------------------------------------------------
pars <- data.frame(
  mu = 0,           # Mean of the predictor (normal distribution)
  sd = 1,           # Standard deviation of the predictor 
  p2 = 0.5,         # Probability of the second predictor (binomial distribution)
  beta_12 = 0.25,   # Effect of the predictor on transition 1 (starting state) to 2 (treatment)
  beta_13 = -0.25,  # Effect of the predictor on transition 1 (starting state) to 3 (recovery)
  beta_23 = -0.15,  # Effect of the predictor on transition 2 (treatment) to 3 (recovery)
  gamma_12 = 0.4,   # Effect of the predictor 2 on transition 1 (starting state) to 2 (treatment)
  gamma_13 = 0.3,   # Effect of the predictor 2 on transition 1 (starting state) to 3 (recovery)
  gamma_23 = 0.1,   # Effect of the predictor 2 on transition 2 (treatment) to 3 (recovery)
  alpha = -0.25,    # Effect of treatment timing on transition 2 (treatment) to 3 (recovery)
  k12 = 0.4,        # Baseline hazard transition 1 (starting state) to 2 (treatment)
  k13 = 0.4,        # Baseline hazard transition 1 (starting state) to 3 (recovery)
  k23 = 0.8         # Baseline hazard transition 2 (treatment) to 3 (recovery)
)

obs_win <- 1.5      # Observation window in yrs
tstrat <- c("0m"=0, "3m"=3/12, "6m"=6/12, "9m"=9/12, "12m"=1, "Never" = 2) # Treatment delay strategies
tps <- seq(0.03, 1.5, 0.03)   # Points at which we want to estimate the outcome
treat_pts <- seq(0, 12, 3)/12 
distr <- "exponential"  # choose between a exponential or weibull (with shape of 0.8)
pred2 <- FALSE          # set true to add a second predictor with a binomial distribution 

# For the transition matrix we use names motivated by our data application:
# Starting state (state 1) is here EM (expectant management)
# Treatment (state 2) is IUI (intrauterine insemination)
# Recovery (state 3) is Preg (pregnancy)
tmat <- trans.illdeath(names=c("EM", "IUI", "preg")) #transition matrix

n <- 2500             # Number of patients
n_sim <- 200           # number of replicates


### Simulation -----------------------------------------------------------------
### ----------------------------------------------------------------------------
sim_data <- function(n, pars, tstrat, time_points, gen_treat, effect_T23, effect_X12, 
                     treat_pts = NULL, multistate_clock = "clock-reset", 
                     distribution = "exponential", predictor2 = FALSE){
  
  ### 0) Generate data ---------------------------------------------------------
  MS <- generate_dat(n, pars, gen_treat = gen_treat, effect_T23 = effect_T23, 
                     effect_X12 = effect_X12, treat_pts = treat_pts, 
                     distribution = distribution, predictor2 = predictor2)
  
  if(!(multistate_clock %in% c("clock-reset", "clock-forward"))) {
    return("mutistate_clock must be either 'clock-reset' or 'clock-forward'")
  }
  
  ### 1) Mstate: ---------------------------------------------------------------
  ###    Prepare data for fitting a multi-state model
  MS$tPreg <- if_else(MS$IUI == 1, MS$stop2, MS$stop1)
  MS$statPreg <- if_else(MS$IUI == 1, MS$Preg2, MS$Preg1)
  
  covs <- c("f_age", if(predictor2) "bin", "stop1")    #covariates in the model
  MSprep <- msprep(time=c(NA, 'stop1', 'tPreg'), status=c(NA,'IUI', 'statPreg'), data=MS, trans=tmat, keep=covs)
  MSdata <- expand.covs(MSprep, covs, append = TRUE, longnames = FALSE)
  
  ### Select covariates in the model, per transition and treatment timing variable
  c12 <- c("f_age.1", if(predictor2) "bin.1")
  c13 <- c("f_age.2", if(predictor2) "bin.2")
  c23 <- c("f_age.3", if(predictor2) "bin.3")
  timecov <- c("stop1.3")
  
  ### Prepare prediction dataset
  MS_pred <- MS |> select(f_age,stop1) #select the covariates in the dataframe MS
  if(predictor2) MS_pred$bin <- MS$bin
  
  ### 1a) Mstate continuous -----------------------------------------------------
  # Fit the model
  type <- "linear"     
  c_fit <- my_mstate(MSdata, c12, c13, c23, timecov, type = type, times = tstrat, clock = multistate_clock)
  # Extract baseline hazard of transition 13 and 23
  baseline <- my_basehaz(c_fit, centered=FALSE, obs_win) # Baseline hazard up until time horizon
  b13 <- baseline[baseline$strata == 'trans=2',]         # Baseline hazard for transition 1 to 3
  b23 <- baseline[grep('trans=3', baseline$strata),]     # Baseline hazard for transition 2 to 3
  
  # Probabilities of recovery over a fine grid of point until time horizon, 
  # for each treatment strategy ("tstrat") 
  CumDensFun2 <- function(treat_time) {
    CumDensFun(MS_pred, c_fit, treat_time, b23, b13, c(c23, timecov), c13, obs_win, time_points, clock = multistate_clock)
  }
  p_ms <- lapply(tstrat, CumDensFun2)
  if(!all(lengths(p_ms) == length(time_points))) return(p_ms) # this stops the simulation run, in case something is wrong
  
  ### 1b) Mstate strata --------------------------------------------------------
  # Fit the model: for the case where treatment time is generate as a continuous variable the 
  # treatment time variable is divided into categories. The categories can be 
  # customized via "grace", the resulting categories are (t^g - grace; t^g + grace).
  # If some time points before time horizon do not fall under any of the generated categories,
  # some extra categories are created as (t^g + grace; t^g' - grace), where t^g is the
  # delay time of strategy g and t^g' is the "next" delay time of strategy g'.
  type <- "strata"  
  grace <- min(diff(tstrat))/2  #same as mstate strata default
  if(gen_treat == "discrete") grace <- 0.001
  c_fit_st <- my_mstate(MSdata, c12, c13, c23, timecov, type = type, tstrat, grace = grace, clock = multistate_clock)
  base_st <- my_basehaz(c_fit_st, centered=FALSE, obs_win) # Baseline hazard up until time horizon
  b13_st <- base_st[base_st$strata == 'trans=2',]          # Baseline hazard for transition 1 to 3
  b23_st <- base_st[grep('trans=3', base_st$strata),]      # Baseline hazard for transition 2 to 3
  
  # Probabilities of recovery over a fine grid of point until time horizon, 
  # for each treatment strategy ("tstrat")
  CumDensFun2 <- function(treat_time) {
    b23_st <- b23_st[grep(paste("_", eval(treat_time), "y", sep = ""), b23_st$strata),]
    CumDensFun(MS_pred, c_fit_st, treat_time, b23_st, b13_st, c23, c13, obs_win, time_points, clock = multistate_clock)
  }
  
  p_ms_st <- lapply(tstrat, CumDensFun2)
  if(!all(lengths(p_ms_st) == length(time_points))) return(list("p_ms_st", p_ms_st)) # this stops the simulation run, in case something is wrong
  
  ### 2) Clone - censor - reweight methods: variant 1 --------------------------
  MS$id <- 1:n         # Create an ID variable
  ccovs <- c("f_age", if(predictor2) "bin")  # Point out the confounder
  
  ### Clone-censor-reweight ----------------------------------------------------
  # This is a first variant of clone-censor-reweight, where the model for the 
  # weights is based on the time-to-treatment mechanism of the original dataset.
  # This means, essentially, that transition S_12(s|X) is estimated.
  # Probability of remaining uncensored given the time-to-treatment survival function
  # Before time - grace: S_12(time|X)
  # Between time - grace and time + gr: S_12(time - gr |X)
  # After time + grace: S_12(time - gr|X) - S_12(time + gr|X)
  # The grace period is only needed for the case where treatment time is generate 
  # as a continuous variable, and is set to 0.001 for the discrete cases.
  
  dfs <- lapply(tstrat, function(t) cens_weight(MS, ccovs, "time-to-treatment", t, grace, rounding = 2))
  
  # Check for infinite weights -------------------------------------------------
  inf_w <- sum(unlist(lapply(dfs, function(df) is.infinite(df$w)))) # Count number of infinite weights
  # Trim if there are any infinite weights
  if (inf_w!= 0) {
    dfs <- lapply(dfs, function(df) {
      trim <- quantile(df$w[!is.infinite(df$w)], probs = 0.975)
      df$w <- if_else(df$w >= trim, trim, df$w) 
      df
    } )
    warning("Infinite weights detected, trimming to 0.975 percentile of
            non-infinite values.") # Warning message in case of infinite weights
  }
  
  ### Kaplan- Meier ------------------------------------------------------------
  p_cl_km <- lapply(dfs, function(df) 1 - surv_km(df, time_points))
  
  ### Cox ----------------------------------------------------------------------
  # We combine all the created weighted datasets (each representing one different 
  # treatment strategy) into one dataset; one extra column is added to keep track 
  # of the treatment strategy followed in each dataset
  dfs2 <- lapply(1:length(tstrat), function(i) cbind(dfs[[names(tstrat)[i]]], Strategy = tstrat[[i]]))
  df_all <- do.call(rbind, dfs2) 
  # clock-reset after treatment
  tstart <- ave(df_all$tstart, df_all$id, df_all$IUI, FUN = first)
  df_all$tstop <- df_all$tstop - tstart 
  df_all$tstart <- df_all$tstart - tstart
  # Interaction between strategy and IUI
  df_all$Strategy1 <- df_all$Strategy * df_all$IUI 
  # Fit the Cox model
  cox_fit <- coxph(Surv(tstart, tstop, Preg) ~ Strategy1 + strata(IUI), 
                   data = df_all, method="breslow", weights = w)
  # Predictions:
  tt <- seq(0.01, 1.5, 0.01)
  b_ccw <- my_basehaz(cox_fit, centered = FALSE, obs_win, tt)
  cox_ccw_fun <- function(t) {
    bt <- b_ccw$hazard[b_ccw$strata == "IUI=0" & b_ccw$time <= t + 10^-6]
    at <- b_ccw$hazard[b_ccw$strata == "IUI=1" & b_ccw$time <= obs_win + 10^-6- t] 
    at <- at * exp(t*cox_fit$coefficients)
    probs <- 1 - exp(-cumsum(c(bt,at)))
    probs <- probs[round(tt*10^5) %in% round(time_points*10^5)]
  }
  p_cl_cox <- lapply(tstrat, cox_ccw_fun)
  
  ### 3) Clone - censor - reweight methods: variant 2 --------------------------
  
  ### Clone-censor-reweight ----------------------------------------------------
  # This is a second variant of clone-censor-reweight, where the model for the 
  # weights is based on the artificial censoring mechanism of each cloned dataset.
  # This means that each clone has a different model for the weights,
  # instead of one model that is applied differently on each clone.
  # We decided to only report in the main article the results of the first variant,
  # but we present the code for this second variant here, together with our results.
  
  dfs_c <- lapply(tstrat, function(t) cens_weight(MS, ccovs, "time-to-censoring", t, grace, rounding = 2))
  
  # Check for infinite weights -------------------------------------------------
  inf_w_c <- sum(unlist(lapply(dfs_c, function(df) is.infinite(df$w))))
  # Trim if inf_w !=0
  if (inf_w_c!= 0) {
    dfs <- lapply(dfs, function(df) {
      trim <- quantile(df$w[!is.infinite(df$w)], probs = 0.975)
      df$w <- if_else(df$w >= trim, trim, df$w) 
      df
    } )
  }
  
  ### Kaplan- Meier ------------------------------------------------------------
  p_cl_km_c <- lapply(dfs_c, function(df) 1 - surv_km(df, time_points))
  
  ### Cox ----------------------------------------------------------------------
  dfs2_c <- lapply(1:length(tstrat), function(i) cbind(dfs_c[[names(tstrat)[i]]], Strategy = tstrat[[i]]))
  df_all_c <- do.call(rbind, dfs2_c) 
  tstart <- ave(df_all_c$tstart, df_all_c$id, df_all_c$IUI, FUN = first)
  df_all_c$tstop <- df_all_c$tstop - tstart #clock-reset
  df_all_c$tstart <- df_all_c$tstart - tstart
  df_all_c$Strategy1 <- df_all_c$Strategy * df_all_c$IUI # Interaction between strategy and IUI
  cox_fit_c <- coxph(Surv(tstart, tstop, Preg) ~ Strategy1 + strata(IUI), 
                     data = df_all_c, method="breslow", weights = w)
  b_ccw_c <- my_basehaz(cox_fit_c, centered = FALSE, obs_win, tt)
  cox_ccw_fun_c <- function(t) {
    bt <- b_ccw_c$hazard[b_ccw_c$strata == "IUI=0" & b_ccw_c$time <= t + 10^-6]
    at <- b_ccw_c$hazard[b_ccw_c$strata == "IUI=1" & b_ccw_c$time <= obs_win + 10^-6- t] 
    at <- at * exp(t*cox_fit_c$coefficients)
    probs <- 1 - exp(-cumsum(c(bt,at)))
    probs <- probs[round(tt*10^5) %in% round(time_points*10^5)]
  }
  p_cl_cox_c <- lapply(tstrat, cox_ccw_fun_c)
  
  
  ### 4) Results ---------------------------------------------------------------
  l <- length(time_points)
  df <- data.frame(
    Method = rep(c("Mstate", "Mstate strata", "CCR + KM", "CCR + Cox", "CCR + KM Cens", "CCR + Cox Cens"), 
                 each = l),
    time = rep(time_points, 6),
    Never = c(p_ms[["Never"]], p_ms_st[["Never"]], p_cl_km[["Never"]], 
              p_cl_cox[["Never"]], p_cl_km_c[["Never"]], p_cl_cox_c[["Never"]]),
    P0m = c(p_ms[["0m"]], p_ms_st[["0m"]], p_cl_km[["0m"]], 
            p_cl_cox[["0m"]], p_cl_km[["0m"]], p_cl_cox[["0m"]]),
    P3m = c(p_ms[["3m"]], p_ms_st[["3m"]], p_cl_km[["3m"]], 
            p_cl_cox[["3m"]], p_cl_km_c[["3m"]], p_cl_cox_c[["3m"]]),
    P6m = c(p_ms[["6m"]], p_ms_st[["6m"]], p_cl_km[["6m"]], 
            p_cl_cox[["6m"]], p_cl_km_c[["6m"]], p_cl_cox_c[["6m"]]),
    P9m = c(p_ms[["9m"]], p_ms_st[["9m"]], p_cl_km[["9m"]], 
            p_cl_cox[["9m"]], p_cl_km_c[["9m"]], p_cl_cox_c[["9m"]]),
    P12m = c(p_ms[["12m"]], p_ms_st[["12m"]], p_cl_km[["12m"]], 
             p_cl_cox[["12m"]], p_cl_km_c[["12m"]], p_cl_cox_c[["12m"]])
  )
  
  row.names(df) <- NULL
  
  return(df)
}

### Run simulation -------------------------------------------------------------
### ----------------------------------------------------------------------------

### Scenario 1) ----------------------------------------------------------------
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="discrete", effect_T23 = "linear", 
           effect_X12 = "linear", treat_pts, distribution = distr, predictor2 = pred2), 
  simplify = F
)
summary_df <- bind_rows(res, .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win, distribution = distr, predictor2 = pred2)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = "Output/Discrete.rda")
#save(summary_df, truth_df, file = "Output/Discrete_wei.rda") # if distr <- "weibull" 
#save(summary_df, truth_df, file = "Output/Discrete_twopred.rda") # if pred2 <- TRUE 

### Scenario 2) ----------------------------------------------------------------
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="continuous", effect_T23 = "linear", 
           effect_X12 = "linear", treat_pts, distribution = distr, predictor2 = pred2), 
  simplify = F
) #30 warnings due to infinite weights
summary_df <- bind_rows(res, .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win, distribution = distr, predictor2 = pred2)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = "Output/continuous_30warn.rda")
#save(summary_df, truth_df, file = "Output/continuous_29warn_wei.rda") # if distr <- "weibull" 
#save(summary_df, truth_df, file = "Output/continuous_31warn_twopred.rda") # if pred2 <- TRUE 

### Scenario 3) ----------------------------------------------------------------
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="discrete", effect_T23 = "non-linear", 
           effect_X12 = "linear", treat_pts, distribution = distr, predictor2 = pred2), 
  simplify = F
)
summary_df <- bind_rows(res, .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win, effect_T23 = "non-linear", 
                                     distribution = distr, predictor2 = pred2)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = "Output/non-linearT.rda")
#save(summary_df, truth_df, file = "Output/non-linearT_wei.rda") # if distr <- "weibull" 
#save(summary_df, truth_df, file = "Output/non-linearT_twopred.rda") # if pred2 <- TRUE 

### Scenario 4) ----------------------------------------------------------------
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="discrete", effect_T23 = "linear", 
           effect_X12 = "non-linear", treat_pts, distribution = distr, predictor2 = pred2), 
  simplify = F
)
summary_df <- bind_rows(res, .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win, distribution = distr, predictor2 = pred2)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = "Output/non-linearX.rda")
#save(summary_df, truth_df, file = "Output/non-linearX_wei.rda") # if distr <- "weibull" 
#save(summary_df, truth_df, file = "Output/non-linearX_twopred.rda") # if pred2 <- TRUE 


### Plots: ---------------------------------------------------------------------
### ----------------------------------------------------------------------------

# Load one scenario and then generate the respective plot
# Scenarios from the main paper
load("Output/Discrete.rda")
#load("Output/continuous_30warn.rda")
#load("Output/non-linearT.rda")
#load("Output/non-linearX.rda")
# Scenarios from the main paper, but with a weibull baseline, rather than exponential
#load("Output/Discrete_wei.rda")
#load("Output/continuous_29warn_wei.rda")
#load("Output/non-linearT_wei.rda")
#load("Output/non-linearX_wei.rda")
# Scenarios from the main paper, but with an extra (binomial) covariate
#load("Output/Discrete_twopred.rda")
#load("Output/continuous_31warn_twopred.rda")
#load("Output/non-linearT_twopred.rda")
#load("Output/non-linearX_twopred.rda")

# Add time 0 to summary
grid <- expand.grid(rep_id = 1:n_sim, Method = unique(summary_df$Method), time = 0)
summary_df2 <- merge(summary_df, grid, by = c("rep_id", "Method", "time"), all = T)
summary_df2[is.na(summary_df2)] <- 0
# Add time 0 to truth
grid <- expand.grid(Strategy = unique(truth_df$Strategy), time = 0)
truth_df2 <- merge(truth_df, grid, by = c("Strategy", "time"), all = T)
truth_df2[is.na(truth_df2)] <- 0

# mean values of the simulation estimates
summ1 <- summary_df2 |>
  select(-c(rep_id)) |>
  group_by(Method, time) |>
  summarise_all(~ mean(.x)) |>
  ungroup() |>
  pivot_longer(cols = -c(Method, time), names_to = "Strategy", values_to = "Mean") 

summ2 <- summary_df2 |>
  select(-c(rep_id)) |>
  group_by(Method, time) |>
  summarise_all(~ quantile(.x, probs = 0.05)) |>
  ungroup() |>
  pivot_longer(cols = -c(Method, time), names_to = "Strategy", values_to = "lowerCI")

summ3 <- summary_df2 |>
  select(-c(rep_id)) |>
  group_by(Method, time) |>
  summarise_all(~ quantile(.x, probs = 0.95)) |>
  ungroup() |>
  pivot_longer(cols = -c(Method, time), names_to = "Strategy", values_to = "upperCI")

summ <- left_join(summ1, summ2, by = c("Method", "time", "Strategy")) |>
  left_join(summ3, by = c("Method", "time", "Strategy")) |>
  left_join(truth_df2, by = c("time", "Strategy"))

summ$Strategy <- factor(
  summ$Strategy,
  levels = c("0m", "3m", "6m", "9m", "12m", "Never"),
  labels = c("0", "0.25", "0.50", "0.75", "1", "Never")
)

summ <- summ %>%
  filter(Method %in% c("Mstate", "Mstate strata", "CCR + KM", "CCR + Cox"))

summ$Method <- factor(
  summ$Method,
  levels = c("Mstate", "CCR + Cox", "Mstate strata", "CCR + KM"), #"CCR + KM Cens", "CCR + Cox Cens"
  labels = c("Multistate continuous", 
             "Clone-censor-reweight continuous",
             "Multistate categorical",
             "Clone-censor-reweight categorical"
             #"Clone-censor-reweight categorical (Cens)", 
             #"Clone-censor-reweight continuous (Cens)"
  )
)

summ <- summ[summ$Strategy %in% c("0", "0.25", "0.50", "0.75", "Never"),]

cols <- c("Never" = "#D41159", "0.25" = "#0C7B0A", "0.50" = "#FFC20A", 
          "0.75" = "#AA11AA", "0" = "#0C7BDC")

ggplot(summ, aes(x = time, color = Strategy)) +
  geom_line(aes(y = Truth), linetype = "solid", linewidth = 1) +
  geom_line(aes(y = Mean), linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(y = Mean, ymin=lowerCI, ymax=upperCI, fill = Strategy), alpha = 0.3, linetype=0)+
  scale_colour_manual(
    name = "Treat at",
    values = cols,
    aesthetics = c("colour", "fill")
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position="bottom", 
    strip.text.x = element_text(size = 15),
    legend.text = element_text(size=15),
    panel.spacing=unit(2.5,"lines"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,0,-10),
    plot.caption = element_text(hjust = 0, size = 15, vjust = 7)
  ) +
  labs(
    x = "Time", 
    y = "Recovery Probability", 
    fill = "Treat at time", 
    color = "Treat at time",
  ) +
  facet_wrap( ~ Method, ncol = 2) 


### Summary tables: ------------------------------------------------------------
### ----------------------------------------------------------------------------

library(xtable)

tble <- function(summary_df, truth_df) {
  
  Time <- c(1.5)
  
  summaryT <- summary_df|>
    pivot_longer(cols = -c(Method, time, rep_id), names_to = "Strategy", values_to = "Estimate") |>
    left_join(truth_df, by = c("time", "Strategy")) 
  
  summT <- summaryT |>
    filter(time %in% Time) |>
    filter(Method %in% c("Mstate", "Mstate strata", "CCR + KM", "CCR + Cox")) |>
    group_by(Strategy, Method, time) |>
    mutate(avg = mean(Estimate)) |>
    summarise(
      Bias = mean(Estimate - Truth),
      SD = (mean((Estimate - avg)^2)^(1/2)),
      RMSE = (mean((Estimate - Truth)^2)^(1/2))
    ) |>
    ungroup() 
  
  summT$Method <- factor(
    summT$Method,
    levels = c("Mstate", "Mstate strata", "CCR + KM", "CCR + Cox", "CCR + KM Cens", "CCR + Cox Cens"),
    labels = c("Multistate continuous", "Multistate categorical",
               "CCR categorical", 
               "CCR continuous",
               "CCR categorical (Cens)", 
               "CCR continuous (Cens)"
    )
  )
  summT$Strategy <- factor(
    summT$Strategy,
    levels = c("0m", "3m", "6m", "9m", "12m", "Never"),
    labels = c("0 months", "3 months", "6 months", "9 months", "12 year", "Never")
  )
  
  summT |>
    select(-c(SD)) |>
    mutate(
      Bias = as.character(round(Bias,2)),
      RMSE = as.character(round(RMSE, 2))
    ) |>
    pivot_longer(cols = c(Bias, RMSE), names_to = "Measure", values_to = "value") |>
    mutate(Method = paste(Method, Measure, sep = " ")) |>
    select(-c(Measure)) |>
    pivot_wider(names_from = "Method", values_from = "value") 
  
}

# Table scenario 1
load(here::here("Simulation","Discrete.rda"))
summ1 <- tble(summary_df, truth_df) |>  add_column(Scenario = "1", .before = "Strategy")

# Table scenario 2
load(here::here("Simulation","continuous_30warn.rda"))
summ2 <- tble(summary_df, truth_df) |>  add_column(Scenario = "2", .before = "Strategy")

# Table scenario 3
load(here::here("Simulation","non-linearT.rda"))
summ3 <- tble(summary_df, truth_df) |>  add_column(Scenario = "3", .before = "Strategy")

load(here::here("Simulation","non-linearX.rda"))
summ4 <- tble(summary_df, truth_df) |>  add_column(Scenario = "4", .before = "Strategy")

final <- rbind(summ1, summ2, summ3, summ4) |>
  relocate(contains("Multistate"), .after = "Strategy")

print(xtable(final), include.rownames=FALSE)


