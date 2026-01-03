prepare_data <- function() {

  library(tidyverse)
  #library(magrittr)
  #library(caTools)
  
  # Read in processed data ----
  dat_count <- readRDS("data/processed/count_data.rds") 
  
  dat_rlen <- readRDS("data/processed/reach_lengths.rds") |>
    select(reach, length) |>
    drop_na()
  
  dat_daily_temperature <- readRDS("data/processed/daily_air_temperature_era5.rds")
  
  dat_daily_discharge <- readRDS("data/processed/daily_discharge.rds")
  
  # Additional data processing ----
  
  ## Survey conditions (visibility and discharge) ----
  dat_conditions <- read_csv("data/raw/conditions.csv") %>%
    select(-comment) %>%
    left_join(dat_daily_discharge, by = "date") %>%
    rename(discharge = discharge.x, discharge_wsc = discharge.y) |>
    mutate(year = year(date)) |>
    summarize(vis = mean(visibility, na.rm = TRUE),
              dis = mean(discharge_wsc, na.rm = TRUE),
              .by = year) |>
    mutate(zvis = as.numeric(scale(vis)),
           zdis = as.numeric(scale(dis)))
  
  # Survey length (number of lanes x reach length) ---
  lreach <- read_csv("data/raw/lane_reach.csv") |>
    left_join(dat_rlen, by = "reach") |>
    group_by(year, date, survey_day) |>
    summarize(lane_reach = sum(nlanes * length)) |>
    mutate(lane_reach = lane_reach / 1000) |> # convert to km
    ungroup() |>
    select(year, survey_day, lane_reach) |>
    pivot_wider(names_from = year, values_from = lane_reach) |>
    select(-survey_day) |>
    as.matrix()
  
  # Number of survey days ----
  sdays <- read_csv("data/raw/lane_reach.csv") |>
    group_by(year) %>%
    summarize(sdays = max(survey_day)) %>%
    pull()
  
  # Proportion of annual interval when survey was conducted ----
  ptm <- read_csv("data/raw/lane_reach.csv") |>
    select(year, date) |>
    mutate(yday = yday(date),
           pspa = ymd(paste(year, 05, 01, sep = "-")),
           yspa = yday(pspa)) |>
    summarize(date = mean(yday),
              pspa = mean(pspa),
              yspa = mean(yspa), .by = year) |>
    # Add 2013 to the data frame and set survey date arbitrarily to the midpoint
    # between 2013 and 2014 (no effect on model as it's not used anyway for 
    # 2013). Also need to add 2023 peak spawning date to compute the proportion 
    # of elapsed time between when the 2022 survey occurred between peak
    # spawning dates of 2022 and 2023.
    bind_rows(tibble(year = c(2013, 2023), 
                     date = yday(ymd(c("2013-05-01", NA))) + 183, 
              pspa = ymd(c("2013-05-01", "2023-05-01")), 
              yspa = yday(ymd(c("2013-05-01", "2023-05-01"))))) |>
    arrange(year) |>
    mutate(elp = as.numeric(difftime(lead(pspa), pspa, units = "days")), 
           ptm = (date - yspa) / elp) |>
    filter(year != 2023) |>
    pull(ptm)
  
  ## Mean summer air temperature ----
  dat_temp <- dat_daily_temperature |>
    mutate(year = year(date)) |>
    filter(!(month(date) == 6 & day(date) < 15)) |>
    filter(!(month(date) == 9 & day(date) > 15)) |>
    summarize(air_temp = mean(air_temp, na.rm = TRUE), 
              .by = year) |>
    mutate(ztair = as.vector(scale(air_temp)))

  ## Mean spring discharge ----
  dat_spr_disch <- dat_daily_discharge |>
    filter(between(month(date), 4, 6), year(date) >= 1988) |>
    filter(!(month(date) == 4 & day(date) < 15)) |>
    arrange(date)|> 
    mutate(year = year(date)) |>
    group_by(year) |>
    # Use na.rm = TRUE because discharge data for 2000 are missing until May 10.
    summarize(dsch = mean(discharge, na.rm = TRUE)) |>
    mutate(zdsch = as.numeric(scale(dsch)))
  
  ## Sockeye salmon abundance ----
  dat_sockeye <- read_csv("data/raw/sk_nuseds.csv") |>
    filter(between(year, 1979, 2022)) |>
    arrange(year) |>
    mutate(mskt2 = lag(abundance),
           mskt1 = abundance,
           zskt2 = as.numeric(scale(mskt2)),
           zskt1 = as.numeric(scale(mskt1))) |>
   filter(year != 1979)
  
  ## Age and length data ----
  dat_age_len <- read_csv("data/raw/age_length.csv") |>
    pivot_longer(`1`:`9`, names_to = "age", values_to = "length") |>
    drop_na(length) |>
    mutate(length = length / 10, # convert mm to cm 
           size_class = case_when(length < 10 ~ NA,
                                  length >= 10 & length < 30 ~ 1,
                                  length >= 30 & length < 50 ~ 2, 
                                  length >= 50 ~ 3),
           age = as.numeric(age),
           year = year - (age_capture - age),
           t = year - 1979) |>
    left_join(select(dat_temp, year, ztair), by = join_by(year)) |>
    left_join(select(dat_sockeye, year, zskt2, zskt1), by = join_by(year))
  
  C <- select(dat_age_len, sample_id, age, size_class) |> 
    pivot_wider(names_from = age, values_from = size_class) |>
    arrange(sample_id) |>
    select(`1`:`9`) |>
    as.matrix()
  
  t_mat <- select(dat_age_len, sample_id, age, t) |> 
    pivot_wider(names_from = age, values_from = t) |>
    arrange(sample_id) |>
    select(`1`:`9`) |>
    as.matrix()
  
  I <- nrow(C)
  
  fA <- apply(C, 1, function(x) min(which(!is.na(x))))
  lA <- apply(C, 1, function(x) max(which(!is.na(x))))
  
  years <- 1988:2022 - 1987
  
  survey_years <- years[-c(26)] # remove 2013 (26)
  
  mu_M <- apply(dat_count, c(2, 3), max)[1, ]
  
  dat <- list(Y = dat_count,
              mu_M = mu_M,
              vis = dat_conditions$vis,
              dis = dat_conditions$dis,
              zvis = dat_conditions$zvis,
              zdis = dat_conditions$zdis,
              d = lreach,
              tair = dat_temp$air_temp,
              dsch = dat_spr_disch$dsch,
              skt2 = dat_sockeye$mskt2,
              skt1 = dat_sockeye$mskt1,
              ztair = dat_temp$ztair,
              zdsch = dat_spr_disch$zdsch,
              zskt2 = dat_sockeye$zskt2,
              zskt1 = dat_sockeye$zskt1,
              w = ptm,
              C = C)
  
  cnt <- list(J = sdays,
              sC = dim(dat_count)[3],
              T = dim(dat_count)[2] + 1, # add 1 for skipped survey year (2013)
              sT = nrow(dat_temp),
              sy = survey_years,
              t_mat = t_mat,
              I = I,
              fA = fA,
              lA = lA)
                
  return(list(dat = dat, cnt = cnt))
  
}