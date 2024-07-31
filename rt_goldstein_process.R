###
### Rt GOLDSTEIN PROCESS
###

# LOAD PACKAGES

# pre set up

# Steps
# 1 setup packages and stan connection


# 90 day analysis
library(concRt)
library(JuliaCall)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(cowplot)
library(scales)
library(latex2exp)
library(tidybayes)
library(stringr)
library(EpiEstim)

# follow https://github.com/igoldsteinh/concRt
julia <- julia_setup()

# load required julia package 
julia_library("concRt")

# set path
julia <- julia_setup(JULIA_HOME = "/home/dthill/juliaup/bin")

# Get command line argument
args <- commandArgs(trailingOnly = TRUE)

# Check if the county name is provided
if (length(args) == 0) {
  stop("No county name provided.")
}

# Extract county name
county_name <- args[1]

# functions
# function
eirr_function <- function(dataframe, log, E_prior, I_prior, R1_prior, Rt_prior){
  
  colnames(dataframe) <- c("Date", "predictor")
  
  # add time change
  dataframe$week <- floor_date(dataframe$Date, unit = "week")
  dataframe$param_change_times <- seq(1:nrow(dataframe))
  param_change_times <- dataframe %>%
    arrange(desc(Date)) %>%
    filter(!duplicated(week))%>%
    select(param_change_times)
  param_change_times <- rev(param_change_times$param_change_times)
  
  # new observations
  dataframe$count <- seq(1:nrow(dataframe))
  obstimes <- as.numeric(dataframe$count) # numeric, triplicate
  
  # data (note if it needs to be log transformed)
  if(log == "Yes"){
    data <- log(dataframe$predictor)
    
  } else if(log == "No"){
    data <- dataframe$predictor
    
  }
  
  # choose to sample from prior or posterior
  priors_only <- FALSE
  
  # choose number of samples, number of chains, and seed
  n_samples <- 20L
  n_chains <- 4L
  seed <- 1L
  
  start_time <- Sys.time()
  
  posterior_samples_eirr <- fit_eirrc(data, 
                                      obstimes, 
                                      param_change_times, 
                                      priors_only, 
                                      n_samples, 
                                      n_chains, 
                                      seed,
                                      E_init_mean = E_prior,
                                      I_init_mean = I_prior,
                                      R1_init_mean = R1_prior,
                                      rt_init_mean = log(Rt_prior))
  # next steps
  # create dataframes of posterior/posterior predictive draws
  posterior_output_eirr <- generate_eirrc(posterior_samples_eirr,
                                          data,
                                          obstimes, 
                                          param_change_times,
                                          seed = seed,
                                          E_init_mean = E_prior,
                                          I_init_mean = I_prior,
                                          R1_init_mean = R1_prior,
                                          rt_init_mean = log(Rt_prior))
  
  # create quantiles of time-varying parameters
  eirr_quantiles <- make_timevarying_quantiles(posterior_output_eirr[[2]])
  
  eirr_rt_quantiles <- eirr_quantiles %>% dplyr::filter(name == "rt_t_values")
  
  end_time <- Sys.time()
  
  end_time - start_time
  
  # check for extra row*
  
  # seems to be an extra observation, maybe because it had an incomplete week
  eirr_rt_quantiles <- eirr_rt_quantiles %>%
    filter(time <= length(unique(dataframe$week)))
  
  # merge to time series
  week <- c(param_change_times, param_change_times, param_change_times)
  eirr_rt_quantiles$param_change_times <- week
  
  # select the 95% quantile and merge in the weeks
  eirr_rt_95 <- eirr_rt_quantiles %>%
    filter(.width == 0.95)
  
  weeks <- dataframe %>%
    select(week, param_change_times)
  
  eirr_rt_95 <- left_join(eirr_rt_95, weeks, by = c("param_change_times")) %>%
    mutate(se_rt = (.upper - .lower) / 3.02) %>%
    select(week, value, .lower, .upper, se_rt) %>%
    rename(mean_rt = value, 
           ll_95_rt = .lower,
           ul_95_rt = .upper)
  
  return(eirr_rt_95)
}

# rt function from epiestim

rt_function <- function(dataframe, mean_si, std_si, weekly){
  
  if(weekly == "Yes"){
    
    # create weekly time windows from df
    Time <- nrow(dataframe)
    t_start <- seq(2, Time-6, by = 7) # starting at 2 as conditional on the past observations
    t_end <- t_start + 6 # adding 6 to get 7-day windows as bounds included in window
    
    # change name of case data column
    colnames(dataframe)[2] <- "case_data"
    colnames(dataframe)[1] <- "Date"
    
    # fill in missing dates
    dataframe <- dataframe %>%
      complete(Date = seq.Date(min(Date), max(Date), by = "day")
      ) %>%
      mutate(I = na.approx(case_data) # interpolate missing values
      ) %>% 
      dplyr::select(Date, I)
    
    # column names change
    colnames(dataframe) <- c("dates", "I")
    
    # calculate Rt
    r0_test <- estimate_R(dataframe,
                          method = "parametric_si",
                          config = make_config(list(
                            mean_si = mean_si,
                            std_si = std_si,
                            t_start = t_start,
                            t_end = t_end
                          )))
    plot(r0_test, "R") 
    
    # Table of dates used
    r_number_dates <- r0_test$dates
    
    # Table of R_0 estimates
    r_number_estim <- r0_test$R %>% 
      # Add date for end of the (weekly) estimation
      mutate(dates = r_number_dates[t_end]) %>% 
      # Extract mean and SD of R estimate
      select(dates, `Mean(R)`, `Std(R)`, `Quantile.0.05(R)`, `Quantile.0.95(R)`)
    
    # extract data for cleaner ggplots
    plot_df <- full_join(dataframe, r_number_estim, by = "dates")
    
    # rename fields
    plot_df$se_rt <- plot_df$`Std(R)` / sqrt(nrow(plot_df))
    plot_df <- plot_df %>%
      rename(date = dates,
             mean_rt = `Mean(R)`,
             ll_95_rt = `Quantile.0.05(R)`,
             ul_95_rt = `Quantile.0.95(R)`
      )%>%
      select(date,  mean_rt, se_rt, ll_95_rt, ul_95_rt)
    plot_df$week <- floor_date(plot_df$date, unit = "week")
    plot_df <- plot_df %>% filter(!is.na(mean_rt))
    return(plot_df)
  } else if(weekly == "No"){
    
    # change name of case data column
    colnames(dataframe)[2] <- "case_data"
    colnames(dataframe)[1] <- "Date"
    
    # fill in missing dates
    dataframe <- dataframe %>%
      complete(Date = seq.Date(min(Date), max(Date), by = "day")
      ) %>%
      mutate(I = na.approx(case_data) # interpolate missing values
      ) %>% 
      dplyr::select(Date, I)
    
    # column names change
    colnames(dataframe) <- c("dates", "I")
    
    # calculate Rt
    r0_test <- estimate_R(dataframe,
                          method = "parametric_si",
                          config = make_config(list(
                            mean_si = mean_si,
                            std_si = std_si
                          )))
    
    # Table of dates used
    r_number_dates <- r0_test$dates
    
    # Table of R_0 estimates
    r_number_estim <- r0_test$R %>% 
      # Add date for end of the (weekly) estimation
      mutate(dates = r_number_dates[t_end]) %>% 
      # Extract mean and SD of R estimate
      select(dates, `Mean(R)`, `Std(R)`, `Quantile.0.05(R)`, `Quantile.0.95(R)`)
    
    # extract data for cleaner ggplots
    plot_df <- full_join(dataframe, r_number_estim, by = "dates")
    
    # rename fields
    plot_df$se_rt <- plot_df$`Std(R)` / sqrt(nrow(plot_df))
    plot_df <- plot_df %>%
      rename(date = dates,
             mean_rt = `Mean(R)`,
             ll_95_rt = `Quantile.0.05(R)`,
             ul_95_rt = `Quantile.0.95(R)`
      )%>%
      select(date,  mean_rt, se_rt, ll_95_rt, ul_95_rt)
    
    return(plot_df)
  }
  
}

# load county data file and make corrections
dat <- readRDS("/home/dthill/Rt_wastewater/Rt_data_county.rds")
dat <- dat %>% mutate(county = ifelse(county == 'St Lawrence', 'St_Lawrence', county),
                      county = ifelse(county == 'New York', 'New_York', county)) %>%
  # remove allegany, and hamilton
  filter(county != "Allegany") %>%
  filter(county != "Hamilton")


# prep data
ww_data <- dat %>%
  filter(county == county_name) %>%
  select(date, sars2.7avg) %>%
  filter(date > min(date)+days(36))

# priors
prior_case_df <- dat %>%
  filter(county == county_name) %>%
  filter(date <= min(ww_data$date) - days(18)) %>%
  arrange(desc(date)) %>%
  head(18) %>%
  arrange(date)

new_df <- prior_case_df %>%
  tail(11) %>%
  mutate(total_c = sum(cases_new.7avg) * 5 * 1/3,
         total_i = sum(cases_new.7avg) * 5 * 2/3)
E_temp <- head(new_df$total_c, 1)
I_temp <- head(new_df$total_i, 1)
R1_temp <- sum(prior_case_df$cases_new.7avg) * 5

# prior rt
Rt_temp <- rt_function(prior_case_df %>%
                         select(date, cases_new.7avg), mean_si = 4, std_si = 1, weekly = "Yes")
Rt_temp <- mean(Rt_temp$mean_rt)

# run analysis
county_eirr <- eirr_function(dataframe = ww_data, log = "Yes",
                             E_prior = E_temp, I_prior = I_temp, R1_prior = R1_temp,
                             Rt_prior = Rt_temp)


# save results
output_file <- paste0('result/', county_name, "_processed_data.rds")
saveRDS(object = county_eirr, file = output_file)


