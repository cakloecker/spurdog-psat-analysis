## -----------------------------------------------------
#
## Script name: 02_data-cleaning.R
#
## Purpose of script: Clean summary and ML data 
# Removes first 24 h and any suspicious data after that to remove tagging and capture/pop-off effects, Adds additional variables
#
## Author: Antonia Kloecker
#
## Institution: Institute of Marine Research (IMR), Norway
#
## E-Mail: antonia.kloecker@hi.no
#
# Date Created: 2023-09-19
#
## Input 
# -   `archives-minutely-summary.rds`
# -   `archives-hourly-summary.rds`
# -   `archives-daily-summary.rds`
## Output 
# -   `archives-minutely-cleaned.csv`
# -   `archives-hourly-cleaned.csv`
# -   `archives-daily-cleaned.csv`
#
## -----------------------------------------------------

## Load packages ########################################
library(patchwork)
library(tidyverse)

## Read Data ############################################
# First execute `01_archive-preprocessing.R` to create the required summary files.
tags <- list.files("data/spurdog/rawdata/")

## Archives
# read in minutely, hourly and daily summaries for all tags
data <- readr::read_rds(file = paste0("D:/development/SharksOnTheMove/PSAT/data/spurdog/processed_data/archives-full.rds"))
mdata <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-minutely-summary.rds")) 
hdata <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-hourly-summary.rds"))
ddata <- readr::read_rds(file = paste0("data/spurdog/processed_data/archives-daily-summary.rds"))
mldata <- readr::read_rds(file = paste0("data/spurdog/processed_data/ML-all.rds"))

## Add more predictors.
mdata <- mdata %>% mutate(
  sharkID = tagID,
  Month = month(Date),
  cohort = case_when(
    tagID %in% c("17P0023", "17P0749", "17P0796", "17P0864", "17P0866") ~ "2019",
    tagID %in% c("20P2042", "20P2044", "20P2045", "20P2046", "20P2047") ~ "2020",
    tagID %in% c("20P2660", "21P0036", "21P0037", "21P0038", "21P0039") ~ "2021",
    tagID %in% c("21P2088", "21P2202", "21P2203", "21P2209") ~ "2022"),
  season2 = case_when(
    Month %in% c(12,1,2) ~ "Dec-Feb",
    Month %in% c(3,4,5) ~ "Mar-May",
    Month %in% c(6,7,8) ~ "Jun-Aug",
    Month %in% c(9,10,11) ~ "Sep-Nov"))
levels(mdata$sharkID) <- 1:length(tags)

hdata <- hdata %>% mutate(
  sharkID = tagID,
  Month = month(Date),
  cohort = case_when(
    tagID %in% c("17P0023", "17P0749", "17P0796", "17P0864", "17P0866") ~ "2019",
    tagID %in% c("20P2042", "20P2044", "20P2045", "20P2046", "20P2047") ~ "2020",
    tagID %in% c("20P2660", "21P0036", "21P0037", "21P0038", "21P0039") ~ "2021",
    tagID %in% c("21P2088", "21P2202", "21P2203", "21P2209") ~ "2022"),
  season2 = case_when(
    Month %in% c(12,1,2) ~ "Dec-Feb",
    Month %in% c(3,4,5) ~ "Mar-May",
    Month %in% c(6,7,8) ~ "Jun-Aug",
    Month %in% c(9,10,11) ~ "Sep-Nov"))
levels(hdata$sharkID) <- 1:length(tags)

ddata <- ddata %>% mutate(
  sharkID = tagID,
  Month = month(Date),
  cohort = case_when(
    tagID %in% c("17P0023", "17P0749", "17P0796", "17P0864", "17P0866") ~ "2019",
    tagID %in% c("20P2042", "20P2044", "20P2045", "20P2046", "20P2047") ~ "2020",
    tagID %in% c("20P2660", "21P0036", "21P0037", "21P0038", "21P0039") ~ "2021",
    tagID %in% c("21P2088", "21P2202", "21P2203", "21P2209") ~ "2022"),
  season2 = case_when(
    Month %in% c(12,1,2) ~ "Dec-Feb",
    Month %in% c(3,4,5) ~ "Mar-May",
    Month %in% c(6,7,8) ~ "Jun-Aug",
    Month %in% c(9,10,11) ~ "Sep-Nov"))
levels(ddata$sharkID) <- 1:length(tags)


## Quick inspection ############################################

# Overview of deployment
hdata %>% group_by(tagID) %>% summarise(n = n(), date_start = range(Date)[1], date_end = range(Date)[2])

# visualise what will be removed (first 24 h for each tag)
test <- hdata %>% group_by(tagID) %>% slice(1:(24*7))%>% 
  mutate(cleaned = ifelse(row_number()%in% 1:24, T, F)) %>% ungroup()
ggplot(test, aes(x=DateTime, y=Depth_median, colour = cleaned))+
  geom_line()+
  scale_colour_manual(values = c("royalblue3", "firebrick3"))+
  scale_y_reverse()+
  scale_x_datetime(date_breaks = "1 day",date_labels = "%d", expand = c(0,0), name = "Date")+
  facet_wrap(~tagID, scales = "free_x")

# Investigate recapture events.

# Tag 20P2660
p1 <- ggplot(mdata %>% dplyr::filter(tagID == "20P2660") %>% 
               filter(DateTime >= ymd_hms("2022-09-20 16:00:00") & DateTime < ymd_hms("2022-09-21 12:00:00")),
             aes(x = DateTime, y = Depth_median)) + 
  geom_point()+
  scale_y_reverse()
  

p2 <- ggplot(mdata %>% dplyr::filter(tagID == "20P2660") %>% 
               filter(DateTime >= ymd_hms("2022-09-20 16:00:00") & DateTime < ymd_hms("2022-09-21 12:00:00")),
             aes(x = DateTime, y = MA_mean)) + 
  geom_point()

p3 <- ggplot(mdata %>% dplyr::filter(tagID == "20P2660") %>% 
               filter(DateTime >= ymd_hms("2022-09-20 16:00:00") & DateTime < ymd_hms("2022-09-21 12:00:00")),
             aes(x = DateTime, y = FS95_cum)) + 
  geom_point()
p1/p2/p3


# Tag 20P2046
p1 <- ggplot(mdata %>% dplyr::filter(tagID == "20P2046") %>% filter(DateTime >= ymd_hms("2021-04-19 23:00:00") & DateTime < ymd_hms("2021-04-20 12:00:00")),
             aes(x = DateTime, y = Depth_median)) + 
  geom_point()+
  scale_y_reverse()

p2 <- ggplot(mdata %>% dplyr::filter(tagID == "20P2046") %>% filter(DateTime >= ymd_hms("2021-04-19 23:00:00") & DateTime < ymd_hms("2021-04-20 12:00:00")),
             aes(x = DateTime , y = MA_mean)) + 
  geom_point()

p3 <- ggplot(mdata %>% dplyr::filter(tagID == "20P2046") %>% filter(DateTime >= ymd_hms("2021-04-19 23:00:00") & DateTime < ymd_hms("2021-04-20 12:00:00")),
             aes(x = DateTime, y = FS95_cum)) + 
  geom_point()
p1/p2/p3

# Tag 21P2088 - constant depth (death/capture? @45 m, 1 am) - gear would be pulled up within 2 days, bottom depth around 75 m in that area
p1 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2088") %>% filter(DateTime >= ymd_hms("2023-07-21 23:00:00") & DateTime < ymd_hms("2023-07-28 12:00:00")),
             aes(x = DateTime, y = Depth_median)) + 
  geom_point()+
  scale_y_reverse()

p2 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2088") %>% filter(DateTime >= ymd_hms("2023-07-21 23:00:00") & DateTime < ymd_hms("2023-07-28 12:00:00")),
             aes(x = DateTime , y = MA_mean)) + 
  geom_point()

p3 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2088") %>% filter(DateTime >= ymd_hms("2023-07-21 23:00:00") & DateTime < ymd_hms("2023-07-28 12:00:00")),
             aes(x = DateTime, y = FS95_cum)) + 
  geom_point()
p1/p2/p3

# Tag 21P2002- got caught in net 25.10.23 (date of pop-up)
p1 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2202") %>% filter(DateTime >= ymd_hms("2023-10-24 23:00:00") & DateTime < ymd_hms("2023-10-25 12:00:00")),
             aes(x = DateTime, y = Depth_median)) + 
  geom_point()+
  scale_y_reverse()

p2 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2202") %>% filter(DateTime >= ymd_hms("2023-10-24 23:00:00") & DateTime < ymd_hms("2023-10-25 12:00:00")),
             aes(x = DateTime , y = MA_mean)) + 
  geom_point()

p3 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2202") %>% filter(DateTime >= ymd_hms("2023-10-24 23:00:00") & DateTime < ymd_hms("2023-10-25 12:00:00")),
             aes(x = DateTime, y = FS95_cum)) + 
  geom_point()
p1/p2/p3


# Tag 21P2003 - check out of interest
p1 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2203") %>% filter(DateTime >= ymd_hms("2023-10-24 23:00:00") & DateTime < ymd_hms("2023-10-27 12:00:00")),
             aes(x = DateTime, y = Depth_median)) + 
  geom_point()+
  scale_y_reverse()

p2 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2203") %>% filter(DateTime >= ymd_hms("2023-10-24 23:00:00") & DateTime < ymd_hms("2023-10-27 12:00:00")),
             aes(x = DateTime , y = MA_mean)) + 
  geom_point()

p3 <- ggplot(mdata %>% dplyr::filter(tagID == "21P2203") %>% filter(DateTime >= ymd_hms("2023-10-24 23:00:00") & DateTime < ymd_hms("2023-10-27 12:00:00")),
             aes(x = DateTime, y = FS95_cum)) + 
  geom_point()
p1/p2/p3


### composite for all recaptures
t1 <- mdata %>% dplyr::filter(tagID == "20P2046") %>% 
  filter(DateTime >= ymd_hms("2021-04-19 15:00:00") & 
           DateTime < ymd_hms("2021-04-20 09:00:00"))
p1 <- ggplot() + 
  geom_rect(aes(xmin = t1$sunset[1], xmax = t1$sunrise[nrow(t1)],
                ymin = min(t1$Depth_median)*0.95, ymax = max(t1$Depth_median)*1.05),
            fill =  "grey60", alpha = 0.3)+
  geom_rect(aes(xmin = t1$dusk.naut[1], xmax = t1$dawn.naut[nrow(t1)],
                ymin = min(t1$Depth_median)*0.95, ymax = max(t1$Depth_median)*1.05),
            fill =  "grey60", alpha = 0.3)+
  geom_point(data = t1, aes(x = DateTime, y = Depth_median), colour = "#5fb667ff")+
  scale_y_reverse()+
  theme_bw()+
  labs(y = "Depth [m]", x = "", title = "Shark 9")

t2 <- mdata %>% dplyr::filter(tagID == "20P2660") %>% 
  filter(DateTime >= ymd_hms("2022-09-20 15:00:00") & DateTime < ymd_hms("2022-09-21 09:00:00"))
p2 <- ggplot() + 
  geom_rect(aes(xmin = t2$sunset[1], xmax = t2$sunrise[nrow(t2)],
                ymin = min(t2$Depth_median)*0.95, ymax = max(t2$Depth_median)*1.05),
            fill =  "grey60", alpha = 0.3)+
  geom_rect(aes(xmin = t2$dusk.naut[1], xmax = t2$dawn.naut[nrow(t2)],
                ymin = min(t2$Depth_median)*0.95, ymax = max(t2$Depth_median)*1.05),
            fill =  "grey60", alpha = 0.3)+
  geom_point(data = t2, aes(x = DateTime, y = Depth_median), colour = "#feb55dff")+
  scale_y_reverse()+
  theme_bw()+
  labs(y = "Depth [m]", x = "", title = "Shark 11")

t3 <- mdata %>% dplyr::filter(tagID == "21P2202") %>% 
  filter(DateTime >= ymd_hms("2023-10-24 15:00:00") & 
        DateTime < ymd_hms("2023-10-25 09:00:00"))

p3 <- ggplot() + 
  geom_rect(aes(xmin = t3$sunset[1], xmax = t3$sunrise[nrow(t3)],
                ymin = min(t3$Depth_median)*0.95, ymax = max(t3$Depth_median)*1.05),
            fill =  "grey60", alpha = 0.3)+
  geom_rect(aes(xmin = t3$dusk.naut[1], xmax = t3$dawn.naut[nrow(t3)],
                ymin = 0.95, ymax = max(t3$Depth_median)*1.05),
            fill =  "grey60", alpha = 0.3)+
  geom_point(data = t3, aes(x = DateTime, y = Depth_median), colour = "#cd0e5bff")+
  scale_y_reverse(expand = c(0,0))+
  theme_bw()+
  labs(y = "Depth [m]", x = "", title = "Shark 17")

p <- p1/p2/p3
ggsave(p, filename = "plots/spurdog/supplement/fig14_recapture_events.png", width = 8, height = 10)

## Data Cleaning ############################################

# First, let's always remove the first 24 h from the timeseries, so we remove any potential tagging effects.
rm_h_prior <- 24 # define number of hours that should be removed
rm_h_post <- 0 # define number of hours that should be removed

mdata <- mdata %>% group_by(tagID) %>% slice((60*rm_h_prior +1):(n()- 60*rm_h_post)) %>% ungroup()
hdata <- hdata %>% group_by(tagID) %>% slice((rm_h_prior +1):(n()-rm_h_post))%>% ungroup()
ddata <- ddata %>% group_by(tagID) %>% slice((rm_h_prior/24 +1):(n()-rm_h_post/24)) %>% ungroup()

# Reduce duration to pop-off date (based on track indicated detachment dates).

hdata <- hdata %>% 
  filter(! (cohort == "2019" & DateTime >= "2020-05-25"),
         ! (tagID == "20P2044" & DateTime >= "2021-06-26"),
         ! (tagID == "20P2045" & DateTime >= "2021-05-30"),
         ! (tagID == "20P2046" & DateTime >= "2021-04-20"),
         ! (tagID == "20P2047" & DateTime >= "2021-07-15"),
         ! (tagID %in% c("21P0036","21P0038","21P0039") & DateTime >= "2022-10-31"),
         ! (tagID %in% c("21P2202","21P2203","21P2209") & DateTime >= "2023-10-25"),
         ! (tagID == "21P2088" & DateTime >= "2023-07-26"),
         ! (tagID == "20P2660" & DateTime >= "2022-09-20"))
mdata <- mdata %>% 
  filter(! (cohort == "2019" & DateTime >= "2020-05-25"),
         ! (tagID == "20P2044" & DateTime >= "2021-06-26"),
         ! (tagID == "20P2045" & DateTime >= "2021-05-30"),
         ! (tagID == "20P2046" & DateTime >= "2021-04-20"),
         ! (tagID == "20P2047" & DateTime >= "2021-07-15"),
         ! (tagID %in% c("21P0036","21P0038" ,"21P0039")  & DateTime >= "2022-10-31"),
         ! (tagID %in% c("21P2202","21P2203","21P2209") & DateTime >= "2023-10-25"),
         ! (tagID == "21P2088" & DateTime >= "2023-07-26"),
         ! (tagID == "20P2660" & DateTime >= "2022-09-20"))
ddata <- ddata %>% 
  filter(! (cohort == "2019" & Date >= "2020-05-25"),
         ! (tagID == "20P2044" & Date >= "2021-06-26"),
         ! (tagID == "20P2045" & Date >= "2021-05-30"),
         ! (tagID == "20P2046" & Date >= "2021-04-20"),
         ! (tagID == "20P2047" & Date >= "2021-07-15"),
         ! (tagID %in% c("21P0036","21P0038" ,"21P0039")  & Date >= "2022-10-31"),
         ! (tagID %in% c("21P2202","21P2203","21P2209") & Date >= "2023-10-25"),
         ! (tagID == "21P2088" & Date >= "2023-07-26"),
         ! (tagID == "20P2660" & Date >= "2022-09-20")) 

# After re-inspection of tags which released prematurely, we remove the following data:
# - tag 17P0866 after 2020-03-21 given a malfunctioning of the tag (e.g. very deep and shallow depth (above surface)).
# - tag 21P0037 after 2022-01-25 given a malfunctioning of the tag.
# - tag 20P2042 after 2021-07-21 due to surfacing.
hdata <- hdata %>% 
  filter(! (tagID == "17P0866" & DateTime >= "2020-03-22"),
         ! (tagID == "21P0037" & DateTime >= "2022-01-26"),
         ! (tagID == "20P2042" & DateTime >= "2021-07-22")) 
mdata <- mdata %>% 
  filter(! (tagID == "17P0866" & DateTime >= "2020-03-22"),
         ! (tagID == "21P0037" & DateTime >= "2022-01-26"),
         ! (tagID == "20P2042" & DateTime >= "2021-07-22")) 
ddata <- ddata %>% 
  filter(! (tagID == "17P0866" & Date >= "2020-03-22"),
         ! (tagID == "21P0037" & Date >= "2022-01-26"),
         ! (tagID == "20P2042" & Date >= "2021-07-22")) 

# We know that the tag sensors have some uncertainty. After a rigorous check, we see that the Depth values remaining do not go below -1.29 m. Since we know that they cannot be # positive and they look like functioning tags, we can set the Depth_median in these cases to zero.
ddata <- ddata %>% mutate(Depth_median = replace(Depth_median, Depth_median < 0, 0))
hdata <- hdata %>% mutate(Depth_median = replace(Depth_median, Depth_median < 0, 0))
mdata <- mdata %>% mutate(Depth_median = replace(Depth_median, Depth_median < 0, 0))


## all in one for raw data: NOTE Error: cannot allocate vector of size 620.6 Mb
data <- data %>% mutate(
  sharkID = tagID,
  Month = month(Date),
  cohort = case_when(
    tagID %in% c("17P0023", "17P0749", "17P0796", "17P0864", "17P0866") ~ "2019",
    tagID %in% c("20P2042", "20P2044", "20P2045", "20P2046", "20P2047") ~ "2020",
    tagID %in% c("20P2660", "21P0036", "21P0037", "21P0038", "21P0039") ~ "2021",
    tagID %in% c("21P2088", "21P2202", "21P2203", "21P2209") ~ "2022"),
  season2 = case_when(
    Month %in% c(12,1,2) ~ "Dec-Feb",
    Month %in% c(3,4,5) ~ "Mar-May",
    Month %in% c(6,7,8) ~ "Jun-Aug",
    Month %in% c(9,10,11) ~ "Sep-Nov")) %>% 
  group_by(tagID) %>% slice((60*12*rm_h_prior +1):(n()-60*12*rm_h_post)) %>% ungroup() %>%
  filter(! (cohort == "2019" & DateTime >= "2020-05-25"),
         ! (tagID == "20P2044" & DateTime >= "2021-06-26"),
         ! (tagID == "20P2045" & DateTime >= "2021-05-30"),
         ! (tagID == "20P2046" & DateTime >= "2021-04-20"),
         ! (tagID == "20P2047" & DateTime >= "2021-07-15"),
         ! (tagID %in% c("21P0036","21P0038","21P0039") & DateTime >= "2022-10-31"),
         ! (tagID %in% c("21P2202","21P2203","21P2209") & DateTime >= "2023-10-25"),
         ! (tagID == "21P2088" & DateTime >= "2023-07-26"),
         ! (tagID == "20P2660" & DateTime >= "2022-09-20"),
         ! (tagID == "17P0866" & DateTime >= "2020-03-22"),
         ! (tagID == "21P0037" & DateTime >= "2022-01-26"),
         ! (tagID == "20P2042" & DateTime >= "2021-07-22")) %>% 
  mutate(Depth = replace(Depth, Depth < 0, 0))

levels(data$sharkID) <- 1:length(tags)

# ml data
mldata <- mldata %>% mutate(
  Month = month(Date),
  cohort = case_when(
    tagID %in% c("17P0023", "17P0749", "17P0796", "17P0864", "17P0866") ~ "2019",
    tagID %in% c("20P2042", "20P2044", "20P2045", "20P2046", "20P2047") ~ "2020",
    tagID %in% c("20P2660", "21P0036", "21P0037", "21P0038", "21P0039") ~ "2021",
    tagID %in% c("21P2088", "21P2202", "21P2203", "21P2209") ~ "2022"),
  season2 = case_when(
    Month %in% c(12,1,2) ~ "Dec-Feb",
    Month %in% c(3,4,5) ~ "Mar-May",
    Month %in% c(6,7,8) ~ "Jun-Aug",
    Month %in% c(9,10,11) ~ "Sep-Nov")) %>% 
  group_by(tagID) %>% slice((rm_h_prior/12+1):(n()-rm_h_post/12)) %>% ungroup() %>% 
  filter(! (cohort == "2019" & DateTime >= "2020-05-25"),
         ! (tagID == "20P2044" & DateTime >= "2021-06-26"),
         ! (tagID == "20P2045" & DateTime >= "2021-05-30"),
         ! (tagID == "20P2046" & DateTime >= "2021-04-20"),
         ! (tagID == "20P2047" & DateTime >= "2021-07-15"),
         ! (tagID %in% c("21P0036","21P0038","21P0039") & DateTime >= "2022-10-31"),
         ! (tagID %in% c("21P2202","21P2203","21P2209") & DateTime >= "2023-10-25"),
         ! (tagID == "21P2088" & DateTime >= "2023-07-26"),
         ! (tagID == "20P2660" & DateTime >= "2022-09-20"),
         ! (tagID == "17P0866" & DateTime >= "2020-03-22"),
         ! (tagID == "21P0037" & DateTime >= "2022-01-26"),
         ! (tagID == "20P2042" & DateTime >= "2021-07-22"))

# Overview of output files
data %>% group_by(tagID) %>% summarise(n = n(), date_start = range(Date)[1], date_end = range(Date)[2])
ddata %>% group_by(tagID) %>% summarise(n = n(), date_start = range(Date)[1], date_end = range(Date)[2])
hdata %>% group_by(tagID) %>% summarise(n = n(), date_start = range(Date)[1], date_end = range(Date)[2])
## Save cleaned summaries ############################################
readr::write_rds(x = data,"data/spurdog/processed_data/archives-full-cleaned.rds", compress = "gz") 
readr::write_rds(x = mdata,"data/spurdog/processed_data/archives-minutely-cleaned.rds", compress = "gz") 
readr::write_rds(x = hdata, "data/spurdog/processed_data/archives-hourly-cleaned.rds", compress = "gz")
readr::write_rds(x = ddata, "data/spurdog/processed_data/archives-daily-cleaned.rds", compress = "gz")
readr::write_rds(x = mldata, "data/spurdog/processed_data/ML-cleaned.rds")
