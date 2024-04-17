### Loading libraries and functions --------------------------------------------
### ----------------------------------------------------------------------------
library(survival)
library(mstate)
library(tidyverse)
library(splines)
library(matrixStats) # for function rowCumsums

source(here::here("Simulation","Simulation Functions.R"))


### Set up parameters ----------------------------------------------------------
### ----------------------------------------------------------------------------
pars <- data.frame(
  mu = 0,           # Mean of the predictor
  sd = 1,           # Standard deviation of the predictor 
  beta_12 = 0.25,   # 0.5 v_na #0.5 v2, # Effect of the predictor on transition EM to IUI
  beta_13 = -0.25,  # -0.5 v_na #-0.75 v2, # Effect of the predictor on transition EM to pregnancy
  beta_23 = -0.15, # -0.25 v_na #-0.25 v2  # Effect of the predictor on transition IUI to pregnancy
  alpha = -0.25,   #-0.15 v1 # Effect of treatment timing, t_A
  k12 = 0.4,       #0.3 v_na  # Baseline hazard transition EM to IUI
  k13 = 0.4,       #0.3 v_na # Baseline hazard transition EM to pregnancy
  k23 = 0.8         # 0.5 v1-0 # Baseline hazard transition IUI to pregnancy
)

obs_win <- 1.5                      # Observation window in yrs
tstrat <- c("0m"=0, "3m"=3/12, "6m"=6/12, "9m"=9/12, "12m"=1, "Never" = 2) # Treatment timing strategies
tps <- seq(0.03, 1.5, 0.03)           # Points at which I want to estimate the outcome: equally spaced, 2 decimals at most
treat_pts <- seq(0, 12, 3)/12 

tmat <- trans.illdeath(names=c("EM", "IUI", "preg")) #transition matrix
# Status 1 = EM; 2 = IUI; 3 = Pregnant
n <- 2500             # Number of patients
n_sim <- 200           # number of replicates

# gen_treat = "discrete"
# effect_T23 ="linear"
# effect_X12 = "non-linear"
# time_points = tps
# multistate_clock = "clock-reset"

### Simulation -----------------------------------------------------------------
### ----------------------------------------------------------------------------
sim_data <- function(n, pars, tstrat, time_points, gen_treat, effect_T23, effect_X12, 
                     treat_pts = NULL, multistate_clock = "clock-reset"){

  ### 0) Generate data ---------------------------------------------------------
  MS <- generate_dat(n, pars, gen_treat = gen_treat, effect_T23 = effect_T23, 
                     effect_X12 = effect_X12, treat_pts = treat_pts)
  
  if(!(multistate_clock %in% c("clock-reset", "clock-forward"))) {
    return("mutistate_clock must be either 'clock-reset' or 'clock-forward'")
  }

  ### 1) Mstate: ---------------------------------------------------------------
  ###    Prepare data for fitting a multi-state model
  MS$tPreg <- if_else(MS$IUI == 1, MS$stop2, MS$stop1)
  MS$statPreg <- if_else(MS$IUI == 1, MS$Preg2, MS$Preg1)
  
  covs <- c("f_age","stop1")          #covariates in the model
  MSprep <- msprep(time=c(NA, 'stop1', 'tPreg'), status=c(NA,'IUI', 'statPreg'), data=MS, trans=tmat, keep=covs)
  MSdata <- expand.covs(MSprep, covs, append = TRUE, longnames = FALSE)
  
  ### Select covariates in the model, per transition and treatment timing variable
  c12 <- c("f_age.1")
  c13 <- c("f_age.2")
  c23 <- c("f_age.3")
  timecov <- c("stop1.3")
  
  ### Prepare prediction dataset
  MS_pred <- MS |> select(f_age,stop1) #select the covariates in the dataframe MS
  
  ### 1a) Mstate continuous -----------------------------------------------------
  type <- "linear"     
  c_fit <- my_mstate(MSdata, c12, c13, c23, timecov, type = type, times = tstrat, clock = multistate_clock)
  baseline <- my_basehaz(c_fit, centered=FALSE, obs_win) # Baseline hazard up until time horizon
  b13 <- baseline[baseline$strata == 'trans=2',]         # Baseline hazard for transition 1 to 3
  b23 <- baseline[grep('trans=3', baseline$strata),]     # Baseline hazard for transition 2 to 3
  
  # Function we apply to all potential treatment points
  CumDensFun2 <- function(treat_time) {
    CumDensFun(MS_pred, c_fit, treat_time, b23, b13, c(c23, timecov), c13, obs_win, time_points, clock = multistate_clock)
  }

  p_ms <- lapply(tstrat, CumDensFun2)
  if(!all(lengths(p_ms) == length(time_points))) return(p_ms)
  
  ### 1b) Mstate strata --------------------------------------------------------
  type <- "strata"     # TODO: implement the type splines
  grace <- min(diff(tstrat))/2  #same as mstate strata default
  if(gen_treat == "discrete") grace <- 0.001
  c_fit_st <- my_mstate(MSdata, c12, c13, c23, timecov, type = type, tstrat, grace = grace, clock = multistate_clock)
  base_st <- my_basehaz(c_fit_st, centered=FALSE, obs_win) # Baseline hazard up until time horizon
  b13_st <- base_st[base_st$strata == 'trans=2',]          # Baseline hazard for transition 1 to 3
  b23_st <- base_st[grep('trans=3', base_st$strata),]      # Baseline hazard for transition 2 to 3
  # Function we apply to all potential treatment points
  CumDensFun2 <- function(treat_time) {
    b23_st <- b23_st[grep(paste("_", eval(treat_time), "y", sep = ""), b23_st$strata),]
    CumDensFun(MS_pred, c_fit_st, treat_time, b23_st, b13_st, c23, c13, obs_win, time_points, clock = multistate_clock)
  }
  
  p_ms_st <- lapply(tstrat, CumDensFun2)
  if(!all(lengths(p_ms_st) == length(time_points))) return(list("p_ms_st",p_ms_st))

  ### 2) Clone - censor - reweight: one treatment model ------------------------
  MS$id <- 1:n
  ccovs <- c("f_age")
  
  ### Reweighting --------------------------------------------------------------
  dfs <- lapply(tstrat, function(t) cens_weight(MS, ccovs, "time-to-treatment", t, grace, rounding = 2))
  
  # Check for infinite weights -------------------------------------------------
  inf_w <- sum(unlist(lapply(dfs, function(df) is.infinite(df$w))))
  # Trim if inf_w !=0
  if (inf_w!= 0) {
    dfs <- lapply(dfs, function(df) {
      trim <- quantile(df$w[!is.infinite(df$w)], probs = 0.975)
      df$w <- if_else(df$w >= trim, trim, df$w) 
      df
      } )
    warning("Infinite weights detected, trimming to 0.975 percentile of
            non-infinite values.")
  }
  
  ### Kaplan- Meier ------------------------------------------------------------
  p_cl_km <- lapply(dfs, function(df) 1 - surv_km(df, time_points))
  
  ### Cox ----------------------------------------------------------------------
  dfs2 <- lapply(1:length(tstrat), function(i) cbind(dfs[[names(tstrat)[i]]], Strategy = tstrat[[i]]))
  df_all <- do.call(rbind, dfs2) 
  tstart <- ave(df_all$tstart, df_all$id, df_all$IUI, FUN = first)
  df_all$tstop <- df_all$tstop - tstart #clock-reset
  df_all$tstart <- df_all$tstart - tstart
  df_all$Strategy1 <- df_all$Strategy * df_all$IUI # Interaction between strategy and IUI
  cox_fit <- coxph(Surv(tstart, tstop, Preg) ~ Strategy1 + strata(IUI), 
                   data = df_all, method="breslow", weights = w)
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
  
  ### 3) Clone - censor - reweight: multiple censoring models ------------------
  ### Reweighting --------------------------------------------------------------
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
           effect_X12 = "linear", treat_pts), 
  simplify = F
  )
summary_df <- bind_rows(res, .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = here::here("Simulation", "Discrete0.rda"))

### Scenario 2) ----------------------------------------------------------------
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="continuous", effect_T23 = "linear", 
           effect_X12 = "linear", treat_pts), 
  simplify = F
)
summary_df <- bind_rows(res[index], .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = here::here("Simulation","continuous0_new.rda"))

### Scenario 3) ----------------------------------------------------------------

### Simulation 
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="discrete", effect_T23 = "non-linear", 
           effect_X12 = "linear", treat_pts), 
  simplify = F
)
summary_df <- bind_rows(res[index], .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win, effect_T23 = "non-linear")

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = here::here("Simulation","non-linearT0.rda"))

### Scenario 4) ----------------------------------------------------------------
set.seed(123456)
res <- replicate(
  n_sim,
  sim_data(n, pars, tstrat, tps, gen_treat ="discrete", effect_T23 = "linear", 
           effect_X12 = "non-linear", treat_pts), 
  simplify = F
)
summary_df <- bind_rows(res, .id = "rep_id")
colnames(summary_df) <- c("rep_id", "Method", "time", "Never", "0m", "3m", "6m", "9m", "12m")

### Truth
truth_fun2 <- function(t)  truth_fun(pars, tps, t, obs_win)

truth_l <- lapply(tstrat, truth_fun2)
truth_df <- do.call(cbind, truth_l) |>
  as.data.frame() |>
  mutate(time = tps) |>
  pivot_longer(cols = -c(time), names_to = "Strategy", values_to = "Truth") 

save(summary_df, truth_df, file = here::here("Simulation","non-linearX_new.rda"))

### Plots: ---------------------------------------------------------------------
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
             "Clone-censor-reweight + Cox",
             "Multistate strata",
             "Clone-censor-reweight + Kaplan-Meier"
             #"Clone-censor-reweight + Kaplan-Meier (Cens)", 
             #"Clone-censor-reweight + Cox (Cens)"
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
    #legend.justification=c(0.975,0),
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
    #caption = "The solid lines represent the truth. The shaded area represents the [5th-95th] percentile of the estimates."
    ) +
  facet_wrap( ~ Method, ncol = 2) 
  #facet_grid(Scenario ~ Method, switch = "y")


### Tables: --------------------------------------------------------------------
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
    labels = c("Multistate continuous", "Multistate strata",
               "CCR + Kaplan-Meier", 
               "CCR + Cox",
               "CCR + Kaplan-Meier (Cens)", 
               "CCR + Cox (Cens)"
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




