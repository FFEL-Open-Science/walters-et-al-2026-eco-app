library(tidyverse)
library(lubridate)
library(tidyhydat)
library(ecmwfr)
library(ncdf4)
download_era5 <- FALSE

# Discharge and water level data from WSC (1988-2022) ----

dat_flows_wsc <- hy_daily_flows(station_number = "08JB002", 
                                start_date = date("1988-01-01"),
                                end_date = date("2022-12-31")) |>
  select(Date, Value) |>
  rename(date = Date, discharge = Value)

saveRDS(dat_flows_wsc, "data/processed/daily_discharge.rds")

# Water temperature from WSC (2011-2022) ----

# These data were obtained upon request from Water Survey of Canada, which 
# does not apply any quality controls check on water temperature data. The 
# checks conducted below were based on the report by Boyer et al (2017) 
# Quality Assured Water Temperature Data from the Water Survey of Canada.
dat_temp_wsc <- read_csv("data/raw/wsc_water_temperature_2011-2024.csv") |>
  arrange(dttm) |> 
  mutate(dttm = update(dttm, minute = 00)) |>
  filter(month(dttm) %in% 6:9, year(dttm) %in% 2011:2022) |>
  group_by(year(dttm)) |>
  complete(dttm = seq(ymd_hms(paste0(unique(year(dttm)), "-06-01 00:00:00")),
                      ymd_hms(paste0(unique(year(dttm)), "-09-30 23:59:00")),
                      by = "hour")) |>
  ungroup() |>
  mutate(dttm_diff = as.numeric(lead(dttm) - dttm),
         abs_temp_lead = abs(lead(temperature) - temperature),
         abs_temp_lag = abs(lag(temperature) - temperature),
         check_1 = between(temperature, 0, 30),
         check_2 = abs_temp_lead < 2 & abs_temp_lag < 2,
         flag_1 = ifelse(check_1 & check_2, "pass", "fail")) |>
  filter(flag_1 == "pass") |>
  group_by(date(dttm)) |>
  mutate(mean_abs_temp_lead = mean(abs_temp_lead),
         mean_abs_temp_lag = mean(abs_temp_lag),
         check_3 = between(mean_abs_temp_lead, 0.01, 0.667) & 
                   between(mean_abs_temp_lag, 0.01, 0.667), 
         flag_2 = ifelse(check_3, "pass", "fail")) |>
  filter(flag_2 == "pass") |>
  mutate(daily_obs = n(),
         check_4 = daily_obs >= 22,
         flag_3 = ifelse(check_4, "pass", "fail"),
         flag_final = ifelse(check_1 & check_2 & check_3 & check_4, 
                             "pass", "fail")) |>
  rename(year = `year(dttm)`, date = `date(dttm)`) |>
  filter(flag_final == "pass") |>
  ungroup()
  
dat_temp_wsc <- summarize(dat_temp_wsc, 
                          water_temp = mean(temperature),
                          .by = date)
  
saveRDS(dat_temp_wsc, "data/processed/daily_water_temperature.rds")

# Air temperature from ERA-5 (1980-2022) ----

if (download_era5 == TRUE) {
  wf_set_key(user = "XXXX", # enter user
             key = "XXXX", # enter key
             service = "cds")
  
  year_start <- c(1980, 1995, 2009)
  year_end <- c(1994, 2008, 2022)
  
  for (i in 1:3) {
    request <- list(
      dataset_short_name = "reanalysis-era5-single-levels",
      product_type = "reanalysis",
      variable = "2m_temperature",
      year = as.character(seq(year_start[i], year_end[i], 1)),
      month = c("06", "07", "08", "09"),
      day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", 
              "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", 
              "23", "24", "25", "26", "27", "28", "29", "30", "31"),
      time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", 
               "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", 
               "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", 
               "21:00", "22:00", "23:00"),
      area = c(54.1, -125, 54, -124.9), # Area including the Stellako
      format = "netcdf",
      target = paste0("era5_air_temperature_", year_start[i], "-", year_end[i], ".nc")
    )
    
    wf_request(user = "321465", request = request, transfer = TRUE,
               path = "data/raw", time_out = 5 * 3600)
  }
}

## Process ERA5 data ----

year_start <- c(1980, 1995, 2009)
year_end <- c(1994, 2008, 2022)

nc_1 <- nc_open("data/raw/era5_air_temperature_1980-1994.nc")
nc_2 <- nc_open("data/raw/era5_air_temperature_1995-2008.nc")
nc_3 <- nc_open("data/raw/era5_air_temperature_2009-2022.nc")

t2m <- c(as.numeric(ncvar_get(nc_1, "t2m")), 
         as.numeric(ncvar_get(nc_2, "t2m")),
         as.numeric(ncvar_get(nc_3, "t2m")))

year_range <- seq(min(year_start), max(year_end), 1)

ny <- length(year_range)

dat_air_temp_era5 <- tibble(
             year = rep(year_range, each = 2928),
             month = rep(c(rep(06, 30 * 24), rep(07, 31 * 24), 
                           rep(08, 31 * 24), rep(09, 30 * 24)), ny),
             day = rep(c(rep(1:30, each = 24), rep(1:31, each = 24), 
                         rep(1:31, each = 24), rep(1:30, each = 24)), ny),
             hour = rep(c(rep(0:23, 30), rep(0:23, 31), 
                          rep(0:23, 31), rep(0:23, 30)), ny)
             ) |>
  mutate(air_temp = t2m - 273.15, # convert to Celsius
         date = ymd(paste0(year, "-", month, "-", day))) |>
  summarize(air_temp = mean(air_temp), .by = date)

saveRDS(dat_air_temp_era5, "data/processed/daily_air_temperature_era5.rds")
