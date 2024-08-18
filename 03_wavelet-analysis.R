## -----------------------------------------------------
##
## Script name: 03_wavelet-analysis.R
##
## Purpose of script: Perform a wavelet analysis on hourly spurdog mean depth and add variable DVM
##
## Author: Antonia Kloecker
##
## Institution: Institute of Marine Research (IMR), Norway
##
## E-Mail: antonia.kloecker@hi.no
##
## Date Created: 2023-09-19
#
## Input 
# -   `archives-hourly-cleaned.csv`
## Output 
# -   `archive_hourly_nested_wavelet.rds`
# -   `archive_hourly_24hwavelet.rds`
## -----------------------------------------------------


## Load packages ########################################
library(tidyverse)
library(WaveletComp)

## Read Data ############################################
hdata <- readr::read_rds("data/spurdog/processed_data/archives-hourly-cleaned.rds")

## Wavelet Analysis #####################################

#depth 1/24
###################
# nest wavelet analysis in tibble
wt24 <- hdata %>% rename(date = DateTime) %>% nest(.by = tagID) %>%
  mutate(wavelet = map(.x = data, ~WaveletComp::analyze.wavelet(my.data = .x,
                                                                my.series = "Depth_median",
                                                                loess.span = 0, # no detrending
                                                                dt = 1/24, # every hour
                                                                # here we choose 1/12 so we have an hour around the 24h and respective signf. levels
                                                                dj = 1/12, #  "resolution" (no of suboctaves) - should correspond with col n.levels
                                                                lowerPeriod = 1/4, # 6 h,
                                                                upperPeriod = 2^7, # 128 d
                                                                # null model is AR1 with p = 0.5
                                                                method = "AR",
                                                                params = list(p = 0.5),
                                                                #make.pval = F,
                                                                n.sim = 1000))) # choose bigger value min. 100 better 1000

readr::write_rds(x = wt24, file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_depth_1_24.rds", compress = "gz")
wt24_red <- wt24 %>% select(-data)
readr::write_rds(x = wt24_red , file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_depth_1_24_wodata.rds", compress = "gz")


# wt24 <-  readr::read_rds(file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_depth_1_24.rds")


# for comparative analysis on same time period
###################
hdata_comp <-  hdata %>% mutate(comp = case_when(sharkID %in% 1:4 & Date >= "2019-12-14" & Date <= "2020-05-23" ~ T,
                                                 sharkID %in% c(6:8,10) & Date >= "2020-12-14" & Date <= "2021-05-24" ~ T,
                                                 sharkID %in% c(11,12,14,15) & Date >= "2021-12-14" & Date <= "2022-05-24" ~ T,
                                                 sharkID %in% 16:19 & Date >= "2022-12-14" & Date <= "2023-05-24" ~ T,
                                                 .default = F))

#hdata_comp %>% filter(comp ==T) %>% group_by(sharkID) %>% summarise(n()) # check that they are same length
wt_comp <- hdata_comp %>% filter(comp ==T) %>% rename(date = DateTime) %>% nest(.by = tagID) %>%
  mutate(wavelet = map(.x = data, ~WaveletComp::analyze.wavelet(my.data = .x,
                                                                my.series = "Depth_median",
                                                                loess.span = 0, # no detrending
                                                                dt = 1/24, #  
                                                                dj = 1/12, #  "resolution" (no of suboctaves) - should correspond with col n.levels
                                                                lowerPeriod = 1/4, # 6 h,
                                                                upperPeriod = 2^7, # 128 d
                                                                # null model is AR1 with p = 0.5
                                                                method = "AR",
                                                                params = list(p = 0.5),
                                                                #make.pval = F,
                                                                n.sim = 100))) # choose bigger value min. 100 better 1000
readr::write_rds(x = wt_comp, file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_comp_1_24.rds", compress = "gz")
#wt_comp <- readr::read_rds(file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_comp_1_24.rds")

#wavelet for fast starts
###########################
wtFS <- hdata %>% rename(date = DateTime) %>% nest(.by = tagID) %>%
  mutate(wavelet = map(.x = data, ~WaveletComp::analyze.wavelet(my.data = .x,
                                                                my.series = "FS95_cum",
                                                                loess.span = 0, # no detrending
                                                                dt = 1/24, # every hour
                                                                # here we choose 1/12 so we have an hour around the 24h and respective signf. levels
                                                                dj = 1/12, #  "resolution" (no of suboctaves) - should correspond with col n.levels
                                                                lowerPeriod = 1/4, # 6 h,
                                                                upperPeriod = 2^7, # 128 d
                                                                # null model is AR1 with p = 0.5
                                                                method = "AR",
                                                                params = list(p = 0.5),
                                                                #make.pval = F,
                                                                n.sim = 1000))) # choose bigger value min. 100 better 1000

readr::write_rds(x = wtFS, file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_FS95cum_1_24.rds", compress = "gz")
#wtFS <- readr::read_rds(file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_FS95cum_1_24.rds")
wtFS_red <- wtFS %>% select(-data)
readr::write_rds(x = wtFS_red , file = "data/spurdog/processed_data/archive_hourly_nested_wavelet_FS95cum_1_24_wodata.rds", compress = "gz")

## Generate scalogramms ####################################

#wt <- readr::read_rds(file = "data/spurdog/processed_data/archive_hourly_nested_wavelet.rds", compress = "gz")


## Generate avg power df to prepare plots  ##############################
get_wtwave <- function(wavelet){return(wavelet$Wave)}
get_wtperiod <- function(wavelet){return(wavelet$Period)}
get_wtpw <- function(wavelet){return(wavelet$Power.avg)}
get_wtpw_sig <- function(wavelet){return(wavelet$Power.avg.pval)}

wtavg <- wt24 %>% 
  #extract average power across periods and respective significance
  mutate(period = map(.x = wavelet, get_wtperiod),
         pw_avg = map(.x = wavelet, get_wtpw),
         pw_avg_sig = map(.x = wavelet, get_wtpw_sig)) %>% 
  # remove wavelet list to avoid replicating them when unnesting
  select(-wavelet) %>% select(-data) %>% 
  unnest(c(period, pw_avg, pw_avg_sig))

readr::write_rds(x = wtavg, file = "data/spurdog/processed_data/archive_hourly_depth_avgwavelet.rds", compress = "gz")
#readr::write_rds(x = wtavg, file = "data/spurdog/processed_data/archive_hourly_depth_FS95cum_avgwavelet.rds", compress = "gz")

### Extract comparive wave data for dendrogramme###########
wtc <- wt_comp %>% mutate(wave = map(.x = wavelet, get_wtwave)) %>% select(-wavelet) %>% select(-data)
readr::write_rds(x = wtc, file = "data/spurdog/processed_data/archive_hourly_depth_compWave.rds", compress = "gz")

## Assess if DVM is present ##############################

#define functions to extract 23-25 h power and significance level
get_24h <- function(wavelet) {
  poi <- which(wavelet$Period > 0.958 & wavelet$Period < 1.042)
  return(wavelet$Power[poi,])
}
get_24h_sig <- function(wavelet) {
  poi <- which(wavelet$Period > 0.958 & wavelet$Period < 1.042)
  return(wavelet$Power.pval[poi,])
}

hdata_wtDVM <- wt24 %>%
  # extract power & resp. significance level at period 23-25 h
  mutate(power24 = map(wavelet, get_24h),
         power24_sig = map(wavelet, get_24h_sig)) %>%
  # remove wavelet list to avoid replicating them when unnesting
  select(-wavelet) %>%
  unnest(c(data, power24, power24_sig)) %>%
  # add logical variable if considered as DVM or not (sig.level <= 0.05)
  mutate(DVM = if_else(power24_sig > 0.05, "nonDVM", "DVM")) %>%
  # reverse the renaming for WaveletComp
  rename(DateTime = date)

readr::write_rds(x = hdata_wtDVM, file = "data/spurdog/processed_data/archive_hourly_DVM_wavelet.rds", compress = "gz")
