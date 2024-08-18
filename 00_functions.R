## -----------------------------------------------------
#
## Script name: 00_functions.R
#
## Purpose of script: Define functions for PSAT spurdog ´04_analysis.Rmd´ script
# for calculating weights for niche plot, plotting time series, wavelet scalogramme etc.
#
## Author: Antonia Kloecker
#
## Institution: Institute of Marine Research (IMR), Norway
#
## E-Mail: antonia.kloecker@hi.no
#
# Date Created: 2023-11-28
#
## -----------------------------------------------------

# Kernel Densities for niche plot ----------------------------------------------

# define function for calculating contours in which x% of points fall
# for plotting contours around a given KD quantile - to use as breaks
density_quantiles <- function(x, y, quantiles) {
  packages <- c("MASS","terra")
  invisible(lapply(packages, library, character.only = TRUE))
  
  dens <- MASS::kde2d(x, y, n = 500)
  df   <- cbind(expand.grid(x = dens$x, y = dens$y), z = c(dens$z))
  r    <- terra::rast(df)
  ind  <- sapply(seq_along(x), function(i) terra::cellFromXY(r, cbind(x[i], y[i])))
  ind  <- ind[order(-r[ind][[1]])]
  vals <- r[ind][[1]]
  ret  <- approx(seq_along(ind)/length(ind), vals, xout = quantiles)$y
  replace(ret, is.na(ret), max(r[]))
}

# obtain weighted densities as required by density_quantiles_w function
kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) 
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}

# for plotting contours around a given KD quantile which is weighted - to use as breaks
density_quantiles_w <- function(x,y,w, quantiles) {
  dens <- kde2d.weighted(x, y, w, n = 150)
  df   <- cbind(expand.grid(x = dens$x, y = dens$y), z = c(dens$z))
  r    <- terra::rast(df)
  ind  <- sapply(seq_along(x), function(i) terra::cellFromXY(r, cbind(x[i], y[i])))
  ind  <- ind[order(-r[ind][[1]])]
  vals <- r[ind][[1]]
  ret  <- approx(seq_along(ind)/length(ind), vals, xout = quantiles)$y
  replace(ret, is.na(ret), max(r[]))
}

# Plotting ----------------------------------------------------------------

#' Plots Depth over time with Temperatures indicated as colours, adds DVM bar on top
#' @param data timeseries data to plot, e.g. hdata, for one individual 
#' @param DVMbar boolean; whether DVM bar should be plotted, default is F
#' @param legend boolean; whether temp legend should be plotted, default is T
#' @param ylims numeric vector with limits for y axis (depth), max depth, min depth, default matches mdata for all spurdog
#' containing covariates Temp_median, Depth_median, DVM with "nonDVM or "DVM", Date and DateTime
plot.DTseries <- function(data, DVMbar = F, legend = T, ylims = c(580,-30)){
  # load required packages & functions
  packages <- c("tidyverse","RColorBrewer","ggnewscale")
  invisible(lapply(packages, library, character.only = TRUE))
    # set up plot with surface line
  p <- ggplot()+
    geom_hline(yintercept=0, linetype='dotted', colour = "grey70")+
  # add depth profile & temp as col
  geom_line(data = data ,
            aes(x = DateTime, y = Depth_median, colour = Temp_median))
  if(legend){
    p <- p + scale_colour_gradientn(colors = c("#013688","#34a2e5ff","#80C586", "#FECC8FFF","#FA815FFF","#bc004cff"), 
                           values = c(0,0.25,0.4,0.6,1),
                           breaks= seq(6,16,2), 
                           limits=c(4,17), 
                           oob = scales::squish, 
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
  }else{
    p <- p + scale_colour_gradientn(colors = c("#013688","#34a2e5ff","#80C586", "#FECC8FFF","#FA815FFF","#bc004cff"), 
                        values = c(0,0.25,0.4,0.6,1),
                        breaks= seq(6,16,2), 
                        limits=c(4,17), 
                        oob = scales::squish, 
                        guide = "none")
    }
  p <- p + scale_y_reverse(expand = c(0,0), limits = ylims, breaks = seq(0,600,100))+
      scale_x_datetime(date_breaks = "1 month", date_labels = "%b",expand = c(0,0))+
  labs(x = "", y = "Depth [m]", colour = "Temp [\u00B0C]")+
  theme_bw()
  
    # add DVM bar
  if(DVMbar){
    # define colour palette of the DVM/nonDVM bar
    myPalette <- colorRampPalette(brewer.pal(3, "Greys"))
    p <- p + ggnewscale::new_scale_colour() +
        geom_line(data = data %>% mutate(DVM = if_else(DVM == "nonDVM", 0,1)), 
            aes(x = DateTime, y = -30, col = DVM), lwd = 8)+
    scale_colour_gradientn(colours = myPalette(2), limits=c(0,1), guide = "none")
  }
  
  return(p)
}


#' Plots Scalogramme from WaveletComp::analyze.wavelet output as a ggplot similar to WaveletComp::wt.image
#' @param wt WaveletComp::analyze.wavelet output
#' @param ts vector; timeseries used for the wavelet analysis in dttm format (including date and time)
#' @param plevel numeric; pvalue that should be used to highlight significant powers with contour
#' @param legened boolean; whether power legend should be plotted, default is T
plot.wavelet <- function(wt, ts, plevel, legend = T){
  
  # load tidyverse library
  invisible(lapply("tidyverse", library, character.only = TRUE))
  
  # create long dataframe from power matrix in WaveletComp output
  wt_df_power <- as.data.frame(wt$Power) %>%
    `colnames<-`(format_ISO8601(ts)) %>%
    mutate(Period = wt$Period) %>%
    pivot_longer(-Period, names_to = "DateTime", values_to = "power") %>%
    mutate(DateTime = lubridate::ymd_hms(DateTime), Period = as.numeric(Period))
  
  # create long dataframe from pval matrix in WaveletComp output
  wt_df_pval <- as.data.frame(wt$Power.pval) %>%
    `colnames<-`(format_ISO8601(ts)) %>%
    mutate(Period = wt$Period) %>%
    pivot_longer(-Period, names_to = "DateTime", values_to = "pval") %>%
    mutate(DateTime = lubridate::ymd_hms(DateTime), Period = as.numeric(Period))
  
  # join both together
  wt_df <- left_join(wt_df_power, wt_df_pval, by = join_by(Period, DateTime))
  
  # Create a df for Cone of influence polygon
  coi_df <- bind_cols(DateTime = lubridate::ymd_hms(c(
    format_ISO8601(ts[1]), # repeat start value
    format_ISO8601(ts[1]), # repeat start value
    format_ISO8601(ts),
    format_ISO8601(ts[length(ts)]), # repeat end value
    format_ISO8601(ts[length(ts)])  # repeat end value
  )), # to match y length with +4 values on the corners
  y = 2^wt$coi.2) # to match log2 y scale
  
  # restrict y limits for plotting reasons
  coi_df$y[which(coi_df$y < 0.25)] <- 0.25
  coi_df$y[which(coi_df$y >128)] <- 132
  
  # make plot
  p <- ggplot()+
    # add "heatmap" with powers over time
    geom_tile(data = wt_df, aes(x = DateTime, y = Period, fill=power))+
    # add significance level as contours
    geom_contour(data = wt_df, aes(x = DateTime, y = Period, z=pval), breaks=c(plevel), colour = "grey30", alpha = 0.7)
    if(legend){
      # set colour scale
     p <- p + scale_fill_gradientn(colors = c("#fef4dcff", "#FDE4A6FF","#FECC8FFF","#FEB37BFF","#FD9A6AFF","#FA815FFF","#F4685CFF","#E85362FF","#c83158ff","#bc004cff","#ab114fff"), 
                           breaks= seq(20,100,20),limits = c(0,85),name = "Power     ",oob = scales::squish, 
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
    }else{
      p <- p + scale_fill_gradientn(colors = c("#fef4dcff", "#FDE4A6FF","#FECC8FFF","#FEB37BFF","#FD9A6AFF","#FA815FFF","#F4685CFF","#E85362FF","#c83158ff","#bc004cff","#ab114fff"), 
                           breaks= seq(20,100,20),limits = c(0,85),name = "Power     ",oob = scales::squish, 
                           guide ="none")
    }
    # log2 transform y scale to match days
    p <- p + scale_y_continuous(trans="log2", breaks = c(0.5,1,7,14,28,56), expand = c(0,0))+
    # set x labels
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b", expand = c(0,0))+
    # add Cone of Influence
    geom_polygon(data = coi_df, aes(x = DateTime, y = y), alpha = 0.7, fill = "white")+
    # define labels and style
    labs(x = "", y = "Period [days]")+
    theme_bw()  
  
  return(p)
}


#' function to add different definitions of night time to Time Series plot for resident species i.e. spurdog, where each day has only one sunset/sunrise
#' @param p ggplot, plot to which night rectangles should be added
#' @param data Dataframe, Time Series data joined with information about night start and end and DateTime col
#' @param from Date, marking the beginning of the period of interest
#' @param to Date, marking the end of the period of interest
#' @param night_start character; name of the col name in data marking beginning of night
#' @param night_end character; name of the col name in data marking end of night
#' @param col character; colour of the night rectangles, default is "grey60"
#' @param alpha integer; transparency of the rectangle, default is 0.3
#' @param ymin_var character or 0; name of the variable to define the minimum limit of the rectangle, default is 0
#' @param ymax_var character or 0; name of the variable to define the maximum y limit of the rectangle, default is 0
#' @param ymmin_scale integer; factor with which the ymax value should be multiplied to expand the y limits of the rectangle, default 0.95
#' @param ymax_scale integer; factor with which the ymax value should be multiplied to expand the y limits of the rectangle, default 1.05
add.Night <- function(p, data, from, to, night_start, night_end, 
                          col = "grey60", alpha = 0.3,
                          ymin_var = 0, ymin_scale = 0.95,
                          ymax_var = 0, ymax_scale = 1.05){
  
  if(data %>%  filter(DateTime >= from & DateTime < to) %>% nrow() == 0)
    stop("There is no data available for the indicated time period.")
  
  if (ymin_var != 0){
    ymin_var <- data %>% 
      dplyr::select(DateTime, !!as.symbol(ymin_var)) %>% 
      filter(DateTime >= from & DateTime < to) %>% 
      summarise_if(is.numeric, min, na.rm= T) %>% 
      pull()*ymin_scale
  }
  
  if (ymax_var != 0){
    ymax_var <- data %>% 
      dplyr::select(DateTime, !!as.symbol(ymax_var)) %>% 
      filter(DateTime >= from & DateTime < to) %>% 
      summarise_if(is.numeric, max, na.rm= T) %>% 
      pull()*ymax_scale
  }    
  
  xmin = data %>% 
    filter(DateTime >= from-1 & DateTime < to+1) %>% 
    distinct(!!as.symbol(night_start)) %>% 
    slice(1:(n()-1)) %>%
    pull()
  
  xmax = data %>% 
    filter(DateTime >= from-1 & DateTime < to+1) %>% 
    distinct(!!as.symbol(night_end)) %>% 
    slice(2:n()) %>% # shifted by a day
    pull() 

  if(all(is.na(xmax))) # in case of polar day or night for entire period
    stop(paste0("There is no available ", night_start," and ", night_end, " time for the period of interest."))
  
  if(any(is.na(xmax))){ # in case of polar day or night in part of the period
    xmax <- xmax[-length(xmax)]
    warning(paste0("For the indicated period there are only part of ", night_start," and ", night_end," times available."))
  }
  
  p <- p + geom_rect(aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin_var,
    ymax = ymax_var),
    fill =  col, alpha = alpha)
  return(p)
}

#' Plots timeseries for depths (lineplot, minutely) coloured with temperature and added activity (barplot, hourly) with second axis for time period of interest, saves plot of interest, nights from sunset to sunrise and from nautical dusk and dawn are added where available
#' @param mdata df/tbl_df; minutely data for depth and temp visualisation in lineplot
#' @param hdata df/tbl_df; hourly data for activity visualisation in barplot
#' @param depth chr; variable name as used in mdata for depth
#' @param depth chr; variable name as used in mdata for depth
#' @param depth chr; variable name as used in mdata for depth
#' @param shark numeric; number of shark (1-19)
#' @param from chr; indicating month and day for start date of period to be displayed, format "MM-DD", year is chosen via cohort covariate
#' @param to chr; indicating month and day for end date of period to be displayed, format "MM-DD", year is chosen via cohort covariate
#' @param season chr; indicating season, used for plot title and file name of saved plot
#' @param legend boolian; whether temperature legend should be ploted, default = T
#' @param filepath chr; filepath to folder in which plot should be saved
#'@description
#'Note: requires a covariate named "cohort" with year of tagging, from and to date are selected in the following year, and "sunset", "sunrise", "naut.dawn", "naut.dusk" as datetimes as provided by Rchival package, if these are not available (NA) no night will be drawn as polygons (based on add.Night function)
#'@example p <- plot.DTAseries(mdata, hdata, shark = 15, from = "06-14", to = "06-21", season = "Summer")
 plot.DTAseries <- function(mdata, hdata,  depth = "Depth_median", temp = "Temp_median",
                           act = "FS95_cum", shark = 1, from = "05-14", to = "05-21", 
                           season = "Summer", legend = T, filepath = "plots/spurdog/supplement/"){
  # load required packages & functions
  packages <- c("tidyverse")
  invisible(lapply(packages, library, character.only = TRUE))
  
  # filter data for shark of interest
  tmp_data <- mdata %>% dplyr::filter(sharkID == shark)
  tmp_data2 <- hdata %>% dplyr::filter(sharkID == shark)
  # define time period of interest
  from <- as.Date(paste0(as.numeric(unique(tmp_data$cohort))+1, "-", from))
  to <- as.Date(paste0(as.numeric(unique(tmp_data$cohort))+1, "-", to))
  
  depth <- "Depth_median" # variable name for depth values
  temp <- "Temp_median" # variable name for temperature values
  act <- "FS95_cum"
  # save temporary data within time period of interest
  tmp_data_toi <- tmp_data %>% filter(DateTime >= from & DateTime < to) 
  tmp_data2_toi <- tmp_data2 %>% filter(DateTime >= from & DateTime < to)
  
  # initialise plot with horizontal line at surface
  p <- ggplot()+geom_hline(yintercept=0, linetype='dotted', colour = "grey70")
  
  # calculate ratio between max values of y axes
  ratio <- max(tmp_data_toi$Depth_median)/max(tmp_data2_toi$FS95_cum)
  # check which y value is bigger - required for appropriate plotting of ydim of the night polygon
  bigger_y <- which.max(c(max(tmp_data_toi$Depth_median),max(tmp_data2_toi$FS95_cum)))
  if(bigger_y == 1){ # depth has max y value
    # check if sunset and sunrise is available (not always the case at high lat)
    if(!any(is.na(tmp_data_toi$sunset)) | !any(is.na(tmp_data_toi$sunrise))){
    p <- add.Night(p = p, data = tmp_data,
                          from = from, to = to,
                          night_start = "sunset",
                          night_end = "sunrise",
                          ymax_var = depth,
                          alpha = 0.2)
    }
    # check if nautical dawn or dusk is available (not always the case at high lat)
    if(!any(is.na(tmp_data_toi$dawn.naut)) | !any(is.na(tmp_data_toi$dusk.naut))){
      p <- add.Night(p = p, 
                       data = tmp_data,
                       from = from, to = to,
                       night_start = "dusk.naut",
                       night_end = "dawn.naut",
                       ymax_var = depth,
                       alpha = 0.2)
    }
    
  }else{ # if act has max y value
    # check if sunset and sunrise is available (not always the case at high lat)
    if(!any(is.na(tmp_data2_toi$sunset)) | !any(is.na(tmp_data2_toi$sunrise))){
    p <- add.Night(p = p, data = tmp_data2,
                          from = from, to = to,
                          night_start = "sunset",
                          night_end = "sunrise",
                          ymax_var = act,
                          ymax_scale = ratio, # scale by ratio to match y axis
                          alpha = 0.2)
    }
    # check if nautical dawn or dusk is available (not always the case at high lat)
    if(!any(is.na(tmp_data2_toi$dawn.naut)) | !any(is.na(tmp_data2_toi$dusk.naut))){
     p <- add.Night(p = p, 
                  data = tmp_data2,
                  from = from, to = to,
                  night_start = "dusk.naut",
                  night_end = "dawn.naut",
                  ymax_var = act,
                  ymax_scale = ratio, # scale factor should be ratio to adjust to y axis
                  alpha = 0.2)
    }
  }
  
  p <- p + 
    # barplot with act
    geom_bar(data = tmp_data2_toi,
                          aes(x = DateTime, y = !!as.symbol(act)*ratio),stat="identity", linewidth = 0.2,  colour = "grey50",fill = "grey70", group = 1)+
    # line plot with depth and temp as col
    geom_line(data = tmp_data_toi,
              aes(x = DateTime, y = !!as.symbol(depth), colour = !!as.symbol(temp)), linewidth = 0.8) +
    scale_x_datetime(date_breaks = "1 day", date_labels = "%d %b",expand = c(0,0))+
    labs(x = "", y = "Depth [m]", colour = "Temp [\u00B0C]", title = season)+
    theme_minimal()+
    scale_y_reverse(sec.axis = sec_axis(~ ./ratio, name = "No. fast starts"), expand = c(0,0))+
    theme_classic()+
    theme(
      axis.title.y = element_text(color = "black"),
      axis.title.y.right = element_text(color = "grey50")
    )
  
  if(legend){
    p <- p + scale_colour_gradientn(colors = c("#013688","#34a2e5ff","#80C586", "#FECC8FFF","#FA815FFF","#bc004cff"), 
                           values = c(0,0.25,0.4,0.6,1),breaks= seq(6,16,2), limits=c(4,17), oob = scales::squish, 
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
  }else{
    p <- p + scale_colour_gradientn(colors = c("#013688","#34a2e5ff","#80C586", "#FECC8FFF","#FA815FFF","#bc004cff"), 
                           values = c(0,0.25,0.4,0.6,1),breaks= seq(6,16,2), limits=c(4,17), oob = scales::squish, 
                           guide = "none")
  }
  ggsave(plot = p, filename = paste0(filepath, "timeseries_", season, "_shark", shark, ".png"), width = 16, height = 12, units = "cm")
  return(p)
}
