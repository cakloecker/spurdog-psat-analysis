#
## Script name: 01_archive-preprocessing.R
#
## Purpose of script: Read and aggregate PSAT archival data to daily, hourly, and 1-minute level
#
## Author: Antonia Kloecker
#
## Institution: Institute of Marine Research (IMR), Norway
#
## E-Mail: antonia.kloecker@hi.no
#
# Date Created: 2023-11-28
#
## Input: 
# -   raw `out-Archive.csv` files
## Output: 
# -   `-archive.rds`
# -   `-archive_solarlunar.rds`
# -   `- archives-daily-summary.rds`
# -   `archives-hourly-summary.rds`
# -   `archives-minutely-summary.rds`
### -----------------------------------------------------

# read packages -----------------------------------------

library(readr)
library(tidyverse)
library(lunar)
library(RchivalTag)

# add solar & lunar params to archive data and calculate summaries --------

tags <- list.files("data/spurdog/rawdata/")


# mixed layer
mldata <- data.frame()
s <- 1
for(t in tags){
  ml_t <- readr::read_csv(paste0("data/spurdog/rawdata/",t, "/out-MixLayer.csv"), col_names = TRUE) %>%
    separate(Date, into = c("Time", "Date"), extra = 'merge', sep = " ") %>% 
    unite("DateTime", c("Date", "Time"), sep = " ", remove = F) %>% 
    mutate(tagID = t,
           sharkID = s,
           Time = lubridate::hms(Time), 
           Date = lubridate::dmy(Date), 
           DateTime = lubridate::dmy_hms(DateTime)) %>% 
    select(tagID, sharkID, DateTime, Date, Time, PerCentMLTime:DepthMax)
  mldata <- rbind(mldata, ml_t)
  s <- s+1
}
readr::write_rds(file = paste0("data/spurdog/processed_data/ML-all.rds"),x = mldata)


# archival data
for(t in tags){
  data <- readr::read_csv(paste0("data/spurdog/rawdata/",t, "/out-Archive.csv"), col_names = TRUE) %>%
    # modify time so it can be read as datetime
    separate(Time, into = c("Time", "Date"), extra = 'merge', sep = " ") %>%
    unite("DateTime", c("Date", "Time"), sep = " ", remove = F) %>%
    mutate(Time = lubridate::hms(Time),
           Date = lubridate::dmy(Date),
           DateTime = lubridate::dmy_hms(DateTime)) %>%
    filter(!if_all(Depth:`Smoothed Light Level`, ~ is.na(.))) %>% # drop rows which only contain NA is all col between Depth and SmoothedLightLevel
    drop_na(DateTime:Temperature) %>% # drops rows where any of the col between DateTime and Temp is NA
    select(!where(~all(. == 0)|all(is.na(.)))) # remove cols which are all NA or zero
  readr::write_rds(file = paste0("data/spurdog/processed_data/",t, "-archive.rds"),x = data)

  #data <- readr::read_rds(file = paste0("data/spurdog/processed_data/",t, "-archive.rds"))

  # add sun and moon dynamics for each day (we can assume the conditions do not change much over one day - increases computation speed)
  data2 <- data %>% distinct(Date, .keep_all = TRUE) %>%
    mutate(Lon = 5.2, Lat = 60.3,# let's assume the sharks stayed close to the tagging location
           datetime = Date + lubridate::hms("12:00:00")) %>% # for consistency, we calculate dawn/dusk also at 12 UTC
    RchivalTag::get_DayTimeLimits() %>%
    select(-datetime) %>%
    mutate(lunar.phase.rad = lunar::lunar.phase(Date), # lunar values refer always to theoretical conditions at 12 UT
           lunar.phase.cat4 = lunar::lunar.phase(Date, name = T),
           lunar.phase.cat8 = lunar::lunar.phase(Date, name = 8),
           lunar.dist = lunar::lunar.distance(Date),
           lunar.illum = lunar::lunar.illumination(Date, shift = 0),
           season = lunar::terrestrial.season(Date),
           dawn.naut = ifelse(is.na(dawn.naut), 100, dawn.naut), # replace NA with 100 to avoid problems with filling
           dawn.ast = ifelse(is.na(dawn.ast), 100, dawn.ast),
           dusk.naut = ifelse(is.na(dusk.naut), 100, dusk.naut),
           dusk.ast = ifelse(is.na(dusk.ast), 100, dusk.ast))

  # # add new covariates to data, copy values for each day and add daynight categories
  data <-  left_join(data, data2) %>%
    fill(Lon, Lat,
         sunrise, sunset, dawn.naut, dawn.ast, dusk.naut, dusk.ast,
          lunar.phase.rad, lunar.phase.cat4, lunar.phase.cat8, lunar.dist,lunar.illum, season) %>%
            # convert 100 back to NA and ensure datetime format
    mutate(dawn.naut = ifelse(dawn.naut == 100, NA, dawn.naut),
         dawn.ast = ifelse(dawn.ast == 100, NA, dawn.ast),
         dusk.naut = ifelse(dusk.naut == 100, NA, dusk.naut),
         dusk.ast = ifelse(dusk.ast == 100, NA, dusk.ast),
         dawn.naut = lubridate::as_datetime(dawn.naut),
         dawn.ast =  lubridate::as_datetime(dawn.ast),
         dusk.naut = lubridate::as_datetime(dusk.naut),
         dusk.ast = lubridate::as_datetime(dusk.ast),
          # add day night variables
         daynight = ifelse(DateTime > sunrise & DateTime < sunset,"day", "night"),
         daynighttwilight.naut = case_when(
           DateTime > dawn.naut & DateTime < sunrise | DateTime > sunset & DateTime < dusk.naut ~ "twilight",
           DateTime > sunrise & DateTime < sunset ~ "day",
           DateTime <= dawn.naut | DateTime >= dusk.naut ~ "night"
         ),
         daynighttwilight.ast = case_when(
           DateTime > dawn.ast & DateTime < sunrise | DateTime > sunset & DateTime < dusk.ast ~ "twilight",
           DateTime > sunrise & DateTime < sunset ~ "day",
           DateTime <= dawn.ast | DateTime >= dusk.ast ~ "night"
         ),
         daynightduskdawn.naut = case_when(
           DateTime > dawn.naut & DateTime < sunrise  ~ "dawn",
           DateTime > sunset & DateTime < dusk.naut ~ "dusk",
           DateTime > sunrise & DateTime < sunset ~ "day",
           DateTime <= dawn.naut | DateTime >= dusk.naut ~ "night"
         ),
         daynightduskdawn.ast = case_when(
           DateTime > dawn.ast & DateTime < sunrise  ~ "dawn",
           DateTime > sunset & DateTime < dusk.ast ~ "dusk",
           DateTime > sunrise & DateTime < sunset ~ "day",
           DateTime <= dawn.ast | DateTime >= dusk.ast ~ "night"
         ),
         tagID = as.factor(t),
         # add activity proxies
         MA = sqrt(Ax^2 + Ay^2 + Az^2),# magnitude of accelleration (on raw accellerometer axes)
         FS99 = MA > quantile(MA, prob = 0.99),  # Fast Starts see Write et al 2021
         FS97 = MA > quantile(MA, prob = 0.97),
         FS95 = MA > quantile(MA, prob = 0.95),
         # add vertical speed
         Depth2 = c(Depth[2:length(Depth)],NA),
         Vertical_speed = (Depth - Depth2)/5, # m/s
         FS_vert = abs(Vertical_speed) > quantile(abs(Vertical_speed), prob = 0.95, na.rm = T)
  ) %>%
  select(-Depth2)

  readr::write_rds(file = paste0("data/spurdog/processed_data/",t, "-archive_solarlunar_MA.rds"),x = data)

  assign(paste0("data_",t), readr::read_rds(file = paste0("data/spurdog/processed_data/",t, "-archive_solarlunar_MA.rds")))

  # minutely data summary
  assign(paste0("m_", t),
      mutate(eval(as.symbol(paste0("data_",t))),
             Minute = lubridate::minute(Time),
             Hour = lubridate::hour(Time)) %>%
      group_by(Date, Hour, Minute, sunrise, sunset, dawn.naut, dawn.ast, dusk.naut, dusk.ast,
               lunar.phase.rad, lunar.phase.cat4, lunar.phase.cat8, lunar.dist,
               lunar.illum, season) %>%
      summarise(Depth_median = median(Depth),
                Vertical_speed_mean = mean(Vertical_speed, na.rm = T), # m/s, accounts for ups and downs -> ameliorating possible tag error
                FS_vert_cum = sum(FS_vert),
                Temp_median = median(Temperature),
                Lightlevel_median = median(`Light Level`, na.rm = T),
                MA_mean = mean(MA),
                FS99_cum = sum(FS99),
                FS97_cum = sum(FS97),
                FS95_cum = sum(FS95),
                # note: does not diff between polar night and day, assumes day
                daynight = ifelse(all(is.na(daynight)), "day", names(which.max(table(daynight)))), # problematic if above polar circle!!
                daynighttwilight.naut = ifelse(all(is.na(daynighttwilight.naut)), "day", names(which.max(table(daynighttwilight.naut)))),
                daynighttwilight.ast = ifelse(all(is.na(daynighttwilight.ast)), "day", names(which.max(table(daynighttwilight.ast)))),
                daynightduskdawn.naut = ifelse(all(is.na(daynightduskdawn.naut)), "day", names(which.max(table(daynightduskdawn.naut)))),
                daynightduskdawn.ast = ifelse(all(is.na(daynightduskdawn.ast)), "day", names(which.max(table(daynightduskdawn.ast))))
      ) %>%
      ungroup() %>%
      mutate(DateTime = lubridate::ymd_hm(paste(Date, Hour, Minute)),
             tagID = as.factor(t))
  )

  # hourly data summary (based on raw data)
  assign(paste0("h_", t),
         mutate(eval(as.symbol(paste0("data_",t))),
                Minute = lubridate::minute(Time),
                Hour = lubridate::hour(Time)) %>%
         group_by(Date, Hour, sunrise, sunset, dawn.naut, dawn.ast, dusk.naut, dusk.ast,
             lunar.phase.rad, lunar.phase.cat4, lunar.phase.cat8, lunar.dist,
             lunar.illum, season) %>%
         summarise(
              Depth_median = median(Depth),
              Vertical_speed_mean = mean(Vertical_speed, na.rm = T) ,# m/s accounts for ups and downs - equal each other out, sign indicates descends and ascends
              Temp_median = median(Temperature),
              Lightlevel_median = median(`Light Level`, na.rm = T),
              MA_mean = mean(MA),
              FS99_cum = sum(FS99),
              FS97_cum = sum(FS97),
              FS95_cum = sum(FS95),
              FS_vert_cum = sum(FS_vert),
              # note: does not diff between polar night and day, assumes day
              daynight = ifelse(all(is.na(daynight)), "day", names(which.max(table(daynight)))), # problematic if above polar circle!!
              daynighttwilight.naut = ifelse(all(is.na(daynighttwilight.naut)), "day", names(which.max(table(daynighttwilight.naut)))),
              daynighttwilight.ast = ifelse(all(is.na(daynighttwilight.ast)), "day", names(which.max(table(daynighttwilight.ast)))),
              daynightduskdawn.naut = ifelse(all(is.na(daynightduskdawn.naut)), "day", names(which.max(table(daynightduskdawn.naut)))),
              daynightduskdawn.ast = ifelse(all(is.na(daynightduskdawn.ast)), "day", names(which.max(table(daynightduskdawn.ast))))
              # this would circumvent the issue of polar night, but takes too long
              # daynight = case_when(
              #   all(is.na(daynight)) & all(season %in% c("Spring", "Summer")) ~ "day", # account for polar day
              #   all(is.na(daynight)) & all(season %in% c("Winter", "Autumn")) ~ "night", # and polar night
              #   !all(is.na(daynight)) ~ names(which.max(table(daynight))),
              #     .default = NA
              # )
      ) %>%
      ungroup() %>%
      mutate(DateTime = lubridate::ymd_h(paste(Date, Hour)),
           tagID = as.factor(t))
    )
  # add vertical absolute speed which is based on minutely means to account for ups and downs created by tag error of +-0.5m
  assign(paste0("h_", t),
         left_join(
           x = eval(as.symbol(paste0("h_", t))),
           y = eval(as.symbol(paste0("m_",t))) %>%
               group_by(Date, Hour) %>%
               summarise(Vertical_speed_abs_mean = mean(abs(Vertical_speed_mean))),
           by = c("Date", "Hour"))
  )

  # daily data summary (based on rawdata)
  assign(paste0("d_", t),eval(as.symbol(paste0("data_",t))) %>%
    group_by(Date, sunrise, sunset, dawn.naut, dawn.ast, dusk.naut, dusk.ast,
             lunar.phase.rad, lunar.phase.cat4, lunar.phase.cat8, lunar.dist,
             lunar.illum, season) %>%
    summarise(
      Depth_median = median(Depth),
      Vertical_speed_mean = mean(Vertical_speed, na.rm = T) ,# m/s accounts for ups and downs - equal each other out, sign indicates descends and ascends
      FS_vert_cum = sum(FS_vert),
      Temp_median = median(Temperature),
      Lightlevel_median = median(`Light Level`, na.rm = T),
      MA_mean = mean(MA),
      FS99_cum = sum(FS99),
      FS97_cum = sum(FS97),
      FS95_cum = sum(FS95)) %>%
    ungroup() %>%
    mutate(tagID = as.factor(t))
    )
  # add vertical absolute speed which is based on minutely means to account for ups and downs created by tag error of +-0.5m
  assign(paste0("d_", t),
         left_join(
           x = eval(as.symbol(paste0("d_", t))),
           y = eval(as.symbol(paste0("m_",t))) %>%
             group_by(Date) %>%
             summarise(Vertical_speed_abs_mean = mean(abs(Vertical_speed_mean))),
           by = "Date")
  )
}

# join and save summaries for all tags ------------------------------------

# # generate vector with element names
d_names <- vector()
h_names <- vector()
m_names <- vector()
names <- vector()
i <- 1
for(t in tags){
  d_names[i] <- paste0("d_", t)
  h_names[i] <- paste0("h_", t)
  m_names[i] <- paste0("m_", t)
  names[i] <- paste0("data_", t)
  i <- i+1
}
# rowbind each data summary across tags
ddata <- data.frame()
hdata <- data.frame()
mdata <- data.frame()
data <- data.frame()

for(i in 1:length(tags)){
  ddata <- rbind(ddata, eval(as.symbol(d_names[i])))
  hdata <- rbind(hdata, eval(as.symbol(h_names[i])))
  mdata <- rbind(mdata, eval(as.symbol(m_names[i])))
  data <- rbind(data, eval(as.symbol(names[i]))) # takes 5-10 min
}

# # save joined summaries
readr::write_rds(file = paste0("data/spurdog/processed_data/archives-minutely-summary.rds"),x = mdata)
readr::write_rds(file = paste0("data/spurdog/processed_data/archives-hourly-summary.rds"),x = hdata)
readr::write_rds(file = paste0("data/spurdog/processed_data/archives-daily-summary.rds"),x = ddata)
readr::write_rds(file = paste0("data/spurdog/processed_data/archives-full.rds"),x = data, compress = "gz")

mdata <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-minutely-summary.rds"))
hdata <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-hourly-summary.rds"))
ddata <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-daily-summary.rds"))
data <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-full.rds"))

unique(hdata$tagID)

