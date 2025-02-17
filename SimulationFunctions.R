### Sequential logistic function -----------------------------------------------
### ----------------------------------------------------------------------------
seq_log <- function(probs, times, horizon = 10) {
  l <- nrow(probs)
  n <- ncol(probs)
  df <- matrix(nrow = l, ncol = n)
  for (i in 1:n)  df[,i] <- rbinom(l, 1, prob = 1-probs[,i])
  df[df ==0] <- NA
  df2 <- t(t(df)*times)
  df2 <- cbind(df2, 10)
  rowMins(df2, na.rm = TRUE)
}

### Draw times from Weibull distribution ---------------------------------------
### ----------------------------------------------------------------------------
rwei <- function(n, a, b) (-log(runif(n))/a)^(1/b)

### Generating data function ---------------------------------------------------
### ----------------------------------------------------------------------------
generate_dat <- function(n,pars, gen_treat = "continuous", effect_T23 = "linear", 
                         effect_X12 = "linear", treat_pts = NULL, 
                         distribution = "exponential", predictor2 = FALSE){
  
  mu <- pars$mu             # Predictor X mean
  sd <- pars$sd             # Predictor X SD
  p2 <- pars$p2
  beta_12 <- pars$beta_12   # Transition linear coefficient of X 1 -> 2
  beta_13 <- pars$beta_13   # Transition linear coefficient of X 1 -> 3
  beta_23 <- pars$beta_23   # Transition linear coefficient of X 2 -> 3
  gamma_12 <- pars$gamma_12 # Transition linear coefficient of X 1 -> 2
  gamma_13 <- pars$gamma_13 # Transition linear coefficient of X 1 -> 3
  gamma_23 <- pars$gamma_23 # Transition linear coefficient of X 2 -> 3
  alpha <- pars$alpha       # Linear coefficient for delay 
  zeta <- pars$zeta         # Linear coefficient of X for censoring
  eta <- pars$eta           # Linear coefficient of predictor 2 for censoring
  k12 <- pars$k12           # Baseline hazard transition 1 -> 2
  k13 <- pars$k13           # Baseline hazard transition 1 -> 3
  k23 <- pars$k23           # Baseline hazard transition 2 -> 3
  kcens <- pars$kcens       # Baseline hazard for censoring
  
  f_age <- rnorm(n, mu, sd)  # Draw n samples from a standard normal distribution
  bin <- 0
  if(predictor2) bin <- rbinom(n, 1, prob = p2)
  dat <- data.frame(f_age)   # Create dataframe
  if(predictor2) dat$bin <- bin

  dat <- dat |> 
    mutate(	
      hr_12 = exp(beta_12 * f_age + gamma_12 * bin),  # Individual HR for EM to IUI 
      hr_13 = exp(beta_13 * f_age + gamma_13 * bin),  # Individual HR for EM to pregnancy
      hr_cens = exp(zeta * f_age + eta * bin),        # Individual HR for censoring
      lambda_12 = hr_12 * k12,                        # Individual hazard from EM to IUI
      lambda_13 = hr_13 * k13,                        # Individual from EM to pregnancy
      lambda_cens = hr_cens * kcens,                  # Individual from EM to pregnancy
      time_12 = rexp(n,lambda_12),                    # Time to IUI drawn from an exp distribution
      time_13 = rexp(n, lambda_13),                   # Time to pregnancy drawn from an exp distribution
      time_cens = rexp(n, lambda_cens))               # Time to censoring drawn from an exp distribution
  
  if(distribution == "weibull") {
    dat$time_12 = rwei(n,dat$lambda_12, 0.8)
    dat$time_13 = rwei(n, dat$lambda_13, 0.8)
    dat$time_cens = rwei(n, dat$lambda_cens, 0.8)
  }
  
  if (gen_treat == "discrete" & effect_X12 == "linear") {
    gr <- min(diff(treat_pts))/2  # To discretize time, equivalent to rounding 
    treat_pts <- treat_pts[order(treat_pts)]
    l <- length(treat_pts)
    probs <- exp(-dat$lambda_12 %*% t(treat_pts + gr))/exp(-dat$lambda_12 %*% t(c(0,(treat_pts +gr)[-l])))
    if(distribution == "weibull") {
      probs <- exp(-dat$lambda_12 %*% t((treat_pts + gr)^0.8))/exp(-dat$lambda_12 %*% t(c(0,(treat_pts +gr)[-l])^0.8))
    }
    dat$time_12 <- seq_log(probs, treat_pts)
  }
  
  if (gen_treat == "discrete" & effect_X12 == "non-linear") {
    gr <- min(diff(treat_pts))/2  # To discretize time, equivalent to rounding 
    treat_pts <- treat_pts[order(treat_pts)]
    l <- length(treat_pts)
    mtx <- matrix(rep(treat_pts + gr,each = n), nrow = n)
    hr <- exp(8*3*beta_12*f_age*(mtx-1/2)^2 - 4*beta_12*f_age + gamma_12 * bin)
    if(distribution == "weibull")  hr <- hr * mtx^-0.2
    dmtx <- matrix(rep(diff(c(0,treat_pts + gr)),each = n), nrow = n)
    probs <- 1-k12*dmtx*hr
    # probs is actually a hazard rate here and may exceed 1 in very extreme cases
    probs[probs<0] <- 0 
    dat$time_12 <- seq_log(probs, treat_pts)
  }
  
  if (effect_X12 == "non-linear" & gen_treat == "continuous") {
    return ("The combination effect_X12 == non-linear & gen_treat == continuous
            was not implemented in this simulation study")
  }
  
  dat <- dat |> 
    mutate(
      hr_23 = exp(beta_23*f_age + gamma_23*bin + alpha*time_12), # Individual HR for IUI to pregnancy 
      lambda_23 = hr_23 * k23,                                  # Individual hazard from IUI to pregnancy
      time_23 = rexp(n, lambda_23)                              # Time from IUI to pregnancy 
      )
  if(distribution == "weibull")  dat$time_23 <- rwei(n, dat$lambda_23, 0.8)  
  
  if (effect_T23 == "non-linear") {
    shape_int <- abs(2*alpha)
    cuts <- c(seq(0, obs_win, 0.1), Inf)
    int <- cut(dat$time_12, cuts, labels = F, right = F)
    t_int <- 1:(length(cuts)-1)
    l <- shape_int/(length(t_int)-1)
    l <- floor(l*10^5)/10^5
    dat$int <- int
    shape_df <- data.frame(
      a = seq(1+alpha, 1-alpha, sign(-alpha) * l),
      int = t_int
    )
    dat <- left_join(dat, shape_df, by = "int")
    dat$lambda_23 <- exp(beta_23 * dat$f_age + gamma_23 * bin) * k23
    dat$time_23 <- rwei(n, dat$lambda_23, dat$a)
    dat <- select(dat, -c(a, int))
  }

  
  dat <- dat |> 
    mutate(
      start1 = rep(0, length(f_age)),               # All patients start from EM set at time = zero
      stop1 = pmin(time_12, time_13),               # Time to event (IUI or natural pregnancy)
      IUI = 1*(time_12 < time_13),                  # IUI: 0 = no, 1 = yes
      Preg1 = 1-IUI,                                # Pregnancy indicator
      ) |>
    mutate(
      IUI = if_else(time_cens < stop1, 0, IUI),     # Censor IUI 
      Preg1 = if_else(time_cens < stop1, 0, Preg1), # Censor event
      stop1 = if_else(time_cens < stop1, time_cens, stop1), # Censor stop1
      start2 = if_else(IUI == 1, time_12, NA_real_),# IUI start time for the treated group
      stop2 = if_else(IUI == 1, time_12 + time_23, NA_real_),  # time to pregnancy for the treated group
      Preg2 = ifelse(IUI == 1, 1, NA_real_)         # Pregnancy for the treated group 
    ) |>
    mutate(
      Preg2 = if_else(time_cens < stop2, 0, Preg2, NA_real_),    # Censored pregnancy for the treated group 
      stop2 = if_else(time_cens < stop2, time_cens, stop2, NA_real_) # Censored time to pregnancy for the treated group
    ) |>
    select(-c(hr_12, hr_23, hr_13, hr_cens, lambda_12, lambda_13, lambda_23, lambda_cens,
              time_12, time_23, time_13, time_cens))
  
  return(dat)
}   

### Function to fit multi-state model ------------------------------------------
### ----------------------------------------------------------------------------
my_mstate <- function(MSdata, c12, c13, c23, timecov, type = "linear", times, grace, clock) {
  
  if (clock == "clock-forward") survt <- "Surv(Tstart, Tstop, status) ~ "
  if (clock == "clock-reset") survt <- "Surv(time, status) ~ "
  
  times <- times[order(times)]
  
  if (type == "linear") {
    fmla <- formula(
      paste(survt, paste(c(c12, c13, c23, timecov), collapse = "+"), "+ strata(trans)")
      )
  }
  
  if (type == "strata") {
    MSdata <- MSdata |> filter(!(trans == 3 & Tstart > obs_win))
    times <- times[times < obs_win]
    if(missing(grace)) grace <- round(min(diff(times))/2, 2) - 0.01
    trt <- MSdata$stop1.3[MSdata$trans == 3]
    for (i in times) {
      trt <- if_else(trt >= i-grace & trt < i+grace, i, trt)
    }
    
    idx <- (!round(10^3*trt) %in% round(10^3*times))
    df_newt <- data.frame(t1 = c(0,times + grace), t2 = c(times - grace, obs_win)) %>%
      filter((t1 < t2) ) %>%
      mutate(tmean = round(rowMeans(.), 2))
    for (i in 1:nrow(df_newt)) {
      trt[idx] <- if_else(trt >= df_newt$t1[i] & trt < df_newt$t2[i], df_newt$tmean[i], trt)[idx]
    }
    
    MSdata$treatt <- 0
    MSdata$treatt[MSdata$trans == 3] <- trt
    MSdata$trans <- case_when(
      MSdata$trans == 1 ~ "trans=1",
      MSdata$trans == 2 ~ "trans=2",
      MSdata$trans == 3 ~ paste("trans=3_", MSdata$treatt, "y", sep = ""),
    )
    
    fmla <- formula(
      paste(survt, paste(c(c12, c13, c23), collapse = "+"), "+ strata(trans)")
    )
  } 
  
  c_fit <- coxph(fmla, data=MSdata, method="breslow")
  
  return(c_fit)
}

### Function to extract baseline up until a fixed time horizon -----------------
### ----------------------------------------------------------------------------

my_basehaz <- function(coxfit, centered = FALSE, obs_win, timepoints) {
  bs <- basehaz(coxfit, centered=centered) |>   # Baseline cumulative hazard for all transitions
    filter(time < obs_win) |>                     # Filter at time of prediction
    rename(cumhaz = hazard)                       # Rename "hazard"to "cumhaz" for clarity
  
  if(is.null(bs$strata)) {
    bs <- bs |>
      mutate(hazard = diff(c(0, cumhaz))) |>
      filter(hazard !=0)  
    
    if(!missing(timepoints)) {
      bs <- bs |>
        select(-hazard) |>
        merge(data.frame(time = timepoints), all = TRUE) |>
        fill(-time) |>
        mutate(cumhaz = if_else(is.na(cumhaz), 0, cumhaz)) |>
        filter(round(time*10^5) %in% round(timepoints*10^5)) |>
        filter(round(10^5*diff(c(0, time)))!=0) |>
        mutate(hazard = diff(c(0, cumhaz)))
    }
  }

  if(!is.null(bs$strata)) {
    bs <- bs |>
      group_by(strata) |>
      mutate(hazard = diff(c(0, cumhaz))) |>        # Compute actual hazard as difference in the hazard
      ungroup() |>
      filter(hazard !=0)     
    if(!missing(timepoints)) {
      grid <- expand.grid(time = timepoints, strata = unique(bs$strata))
      bs <- bs |>
        select(-hazard) |>
        merge(grid, all = TRUE) |>
        arrange(strata, time) |>
        group_by(strata) |>
        fill(-c(time)) |>
        ungroup() |>
        mutate(cumhaz = if_else(is.na(cumhaz), 0, cumhaz)) |>
        filter(round(time*10^5) %in% round(timepoints*10^5)) |>
        filter(round(10^5*diff(c(0, time)))!=0) |>
        group_by(strata) |>
        mutate(hazard = diff(c(0, cumhaz))) |>
        ungroup()
        
    }
  }
  
  return(bs)
}

### Cumulative density function mstate------------------------------------------
### ----------------------------------------------------------------------------

CumDensFun <- function(data, cox_fit, treat_time, base_treated, base_untreated,
                       coeff_treated, coeff_untreated, obs_win, time_points,
                       clock = "clock-reset") {
  if (!(clock %in% c("clock-reset", "clock-forward"))) {
    return ("Choose 'clock' to be 'clock-reset' or 'clock-forward'")
  }
  
  c_fit <- cox_fit
  b23 <- base_treated
  b13 <- base_untreated
  coef23 <- coeff_treated
  coef13 <- coeff_untreated
  
  if (clock == "clock-reset") {
    idx1 <- (b13$time <= treat_time)                 # index pre-treatment times
    idx2 <- (b23$time < obs_win - treat_time)        # index post-treatment times
    time <- c(b13$time[idx1], treat_time + b23$time[idx2])
    if (sum(idx1) == 0) time <- b23$time[idx2]
    bhaz <- c(b13$hazard[idx1], b23$hazard[idx2])    # baseline hazard 
  }
  if (clock == "clock-forward") {
    idx1 <- (b13$time < treat_time)                  # index pre-treatment times
    idx2 <- (b23$time >= treat_time)                 # index post-treatment times
    time <- c(b13$time[idx1], b23$time[idx2])
    bhaz <- c(b13$hazard[idx1], b23$hazard[idx2])    # baseline hazard
  }
  #Patient-specific HRs with and without treatment
  data$stop1 <- treat_time
  HR13 <- exp(rowSums(t(c_fit$coefficients[coef13]*t(data[, sub("[.]2", "", coef13)])))) #HR untreated
  HR23 <- exp(rowSums(t(c_fit$coefficients[coef23]*t(data[, sub("[.]3", "", coef23)])))) #HR treated
  
  # Patient-specific hazards through time given the chosen treatment strategy
  haz <- HR13 %*% t(bhaz*(time < treat_time)) + HR23 %*% t(bhaz*(time >= treat_time))
  
  # Patient-specific cumulative hazard through time given the chosen treatment strategy
  cumhaz <- rowCumsums(haz) 
  Surv <- exp(-cumhaz)    # Survival through time given the chosen treatment strategy
  P <-  Surv * haz
  #EM <- ifelse(treat_time > obs_win, NA_real_, 1 - rowSums(P[, (time < treat_time)]))
  # Record the probabilities at the times of interest
  P_times <- data.frame(t(rowCumsums(P))) |>     
    mutate(time = time) |>
    merge(data.frame(time = time_points), all = TRUE) |>
    fill(-time)
  P_times[is.na(P_times)] <-  0
  P_times <- P_times |> filter(time %in% time_points)
  
  Probs <- rowMeans(select(P_times, -time))
  
  return(Probs)
}

### Create weighted dataset for KM  --------------------------------------------

# The input data needs the following variables: id, tstart1, tstop1, tstart2, 
# tstop2, IUI, Preg1, Preg2 + some predictors

# Weights_type:
# 1) Time-to-treatment = time-to-treatment model fitted based on the entire original data set
# 2) Time-to-censoring = time-to-censoring model fitted based on the specific clone-censor data
#    with the right treatment strategy 

cens_weight <- function(data, pred_vars, weights_type = "time-to-treatment", treat_time, grace, rounding) {
  
  # Data set of consisting of pre-treatment variables 
  MSkm <- select(data, c(id, all_of(pred_vars), tstart = start1, 
                         tstop = stop1, IUI, Preg = Preg1)) %>%
    mutate_at(vars(tstart, tstop), ~ round(.x*10^rounding)) 
  
  # Data set of consisting of post-treatment variables
  MSkm_2 <- select(data, c(id, all_of(pred_vars), tstart = start2, 
                           tstop = stop2, IUI, Preg = Preg2)) %>%
    mutate_at(vars(tstart, tstop), ~ round(.x*10^rounding))
  
  # Censored dataset 
  time <- round(treat_time*10^rounding)
  gr <- round(grace*10^rounding)
  MSkm_i <- MSkm                  
  MSkm_i$sCens <- 0                        # Initiate censoring status
  MSkm_i$tCens <- MSkm_i$tstop             # Initiate censoring time
  idx <- (MSkm_i$IUI == 1)                 # Treated subjects
  idx2 <- (MSkm_i$tstop < time - gr)       # stop1 before supposed treatment time
  idx3 <- (MSkm_i$tstop > time + gr)       # stop1 after supposed treatment time
  MSkm_i$sCens[(idx & idx2)|idx3] <- 1     # censoring those who are treated outside the allowed time period
  MSkm_i$tCens[idx3] <- time + gr          # set censoring time to end of grace period for those 
  MSkm_i$Preg[MSkm_i$sCens==1] <- 0        # pregnancy status = 0 if patient is censored
  
  # Keep post-treatment information on those who are not censored and event-free
  MSkm_2 <- MSkm_2[MSkm_i$sCens==0 & MSkm_i$Preg==0 & MSkm_i$IUI == 1,]
  MSkm_2$tstop[MSkm_2$tstart==MSkm_2$tstop] <- (MSkm_2$tstop+0.1)[MSkm_2$tstart==MSkm_2$tstop]
  
  # Censoring model
  
  if (weights_type == "time-to-censoring") {
    fmla <- formula(paste("Surv(tCens, sCens) ~ ", paste(pred_vars, collapse = "+" )))
    cox_d <- coxph(fmla, data=MSkm_i, method="breslow")
  }
  
  if (weights_type == "time-to-treatment") {
    fmla <- formula(paste("Surv(tstop, IUI) ~ ", paste(pred_vars, collapse = "+" )))
    cox_d <- coxph(fmla, data=MSkm, method="breslow")
  }
  
  # Baseline hazard
  bh_d <- basehaz(cox_d, centered = FALSE)|> rename(cumhaz = hazard) # baseline cumhazard
  bh_d$haz <- diff(c(0,bh_d$cumhaz))                                 # baseline hazard
  bh_d <- bh_d[bh_d$haz != 0,]                                       # remove useless lines
  
  # Weighted dataset
  MSkm_i$tCens[MSkm_i$tCens==0] <- 10^-5
  MSkm_i$tstop<- MSkm_i$tCens
  MSkm_i$IUI<- 0
  MSkm_i <- select(MSkm_i, c(id, all_of(pred_vars), tstart, tstop, IUI, Preg)) |>
    rbind(MSkm_2) |>
    arrange(id, tstop)
  
  if (weights_type == "time-to-censoring")    cuts <- bh_d$time[bh_d$time!=0]
  if (weights_type == "time-to-treatment") {
    cuts <- c(bh_d$time[bh_d$time!=0 & bh_d$time <= time - gr], time + gr)
  }
  
  MSkm_i <- survSplit(data = MSkm_i, cut = cuts, event = "Preg",
                      start = "tstart", end = "tstop")
  MSkm_i$tstop<- round(MSkm_i$tstop)

  MSkm_i <- left_join(MSkm_i,  select(bh_d, c(tstart = time, haz)), by = "tstart") 
  MSkm_i$haz[is.na(MSkm_i$haz)] <- 0
  
  lp <- cox_d$coefficients
  hr <- exp(rowSums(t(t(MSkm_i[, names(lp)]) * lp)))
  MSkm_i$haz <- hr * MSkm_i$haz
  survT <- ave(MSkm_i$haz, MSkm_i$id, FUN = function(a) exp(-cumsum(a)))
  MSkm_i$w <- 1/survT
  if (weights_type == "time-to-treatment") {
    # Probability of remaining uncensored given the time-to-treatment survival function
    # Before time - gr: S(t)
    # Between time - gr and time + gr: S(time - gr)
    # After time + gr: S(time - gr) - S(time + gr)
    beforeT <- (MSkm_i$tstart <= time - gr)
    grwin <- (MSkm_i$tstart > time - gr & MSkm_i$tstart < time + gr)
    afterT <- (MSkm_i$tstart >= time + gr)
    haz2 <- MSkm_i$haz
    haz2[grwin | afterT] <- 0
    survT2 <- ave(haz2, MSkm_i$id, FUN = function(a) exp(-cumsum(a)))
    MSkm_i$w[grwin] <- (1/survT2)[grwin]
    MSkm_i$w[afterT] <- (1/(survT2 - survT))[afterT]
  }
  MSkm_i <- select(MSkm_i, -c(haz)) |>
    mutate(tstop = if_else(tstop == tstart, tstop + 0.1, tstop)) |>
    mutate_at(vars(tstart, tstop), ~ round(.x/10^rounding, rounding + 1))
  
  return(MSkm_i)
} 

surv_km <- function(df1, time_pts1) {
  df_km <-df1[!is.infinite(df1$w) & !is.na(df1$w),]
  if(nrow(df_km) == nrow(df1)) {  #0
    km <- survfit(Surv(tstart, tstop, Preg) ~ 1, data = df_km, weights = w)
    surv <- summary(km, times = time_pts1, extend = TRUE)$surv
  }
  else surv <- rep(NA, length(time_pts1))
  surv
}
# Note that I need the "extend" = TRUE as sometimes I have no subjects left

truth_fun <- function(pars, time_points, treat_t, obs_win, effect_T23 = "linear", 
                      distribution = "exponential", predictor2 = FALSE) {
  vals <- seq(round(qnorm(0.025),3), round(qnorm(0.975),3), by = 0.001)
  time <- rep(time_points, times = length(vals))
  f_age <- rep(vals, each = length(time_points))
  beta_13 <- pars$beta_13
  beta_23 <- pars$beta_23
  gamma_13 <- pars$gamma_13
  gamma_23 <- pars$gamma_23
  alpha <- pars$alpha
  k13 <- pars$k13
  k23 <- pars$k23
  hr13 <- exp(beta_13 * f_age)
  hr23 <- exp(beta_23 * f_age + alpha * treat_t)
  p2 <- pars$p2
  if (predictor2) {
    hr13_2 <- exp(beta_13 * f_age + gamma_13)
    hr23_2 <- exp(beta_23 * f_age + gamma_23 + alpha * treat_t)
  }
  a13 <- a23 <- 1
  
  if (distribution == "weibull") a13 <- a23 <- 0.8
  
  if(effect_T23 == "non-linear") {
    shape_int <- abs(2*alpha)
    cuts <- c(seq(0, obs_win, 0.1), Inf)
    int <- cut(treat_t, cuts, labels = F, right = F)
    t_int <- 1:(length(cuts)-1)
    l <- shape_int/(length(t_int)-1)
    l <- floor(l*10^5)/10^5
    a23 <- seq(1+alpha, 1-alpha, sign(-alpha) * l)[int]
    hr23 <- exp(beta_23 * f_age)
    if (predictor2) hr23_2 <- exp(beta_23 * f_age  + gamma_23)
  }
  
  probs <- (1 - exp(-k13 * time^a13 * hr13)) * (time <= treat_t) +
    ((1 - exp(-k13 * treat_t^a13 * hr13)) +
       (1 - exp(-k23 * abs(obs_win - treat_t + time_points[1] - rev(time))^a23 * hr23)) * 
       (exp(-k13 * treat_t^a13 * hr13))) * (time > treat_t) 
  
  if (predictor2) {
    probs2 <- (1 - exp(-k13 * time^a13 * hr13_2)) * (time <= treat_t) +
      ((1 - exp(-k13 * treat_t^a13 * hr13_2)) +
         (1 - exp(-k23 * abs(obs_win - treat_t + time_points[1] - rev(time))^a23 * hr23_2)) * 
         (exp(-k13 * treat_t^a13 * hr13_2))) * (time > treat_t) 
  }
  
  df <- data.frame(time = time, pr = probs)|>
    group_by(time) |>
    summarise(pr = sum(pr)/n()) |>
    ungroup() 
  pr <- df$pr
  if (predictor2) {
    df2 <- data.frame(time = time, pr = probs2)|>
      group_by(time) |>
      summarise(pr = sum(pr)/n()) |>
      ungroup() 
    pr <- df$pr * (1-p2) + p2 * df2$pr
  }
  
  return(pr)
}









