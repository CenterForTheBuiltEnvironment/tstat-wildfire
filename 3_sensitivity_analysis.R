# Using smart thermostats to reduce indoor exposure to wildfire fine particulate matter (PM2.5)
# Cite: https://doi.org/10.1016/j.indenv.2025.100088
# ----
# TASK: Use a Monte Carlo simulation to find the model parameters that are influencing indoor PM2.5 concentration/exposure
# Code Authors: Federico Dallo, Thomas Parkinson, Carlos Duarte, Chai Yoon Um, Paul Raftery
# ----

#### SETUP ####

# use pacman for package management
library(pacman)

# load packages
pacman::p_load(tidyverse, here, lubridate, # essential
               scales, patchwork, # plots
               randtoolbox, # sobol
               zoo, # timeseries
               sf, # spatial data
               pbapply, # progress bars
               furrr, future.apply # parallel calculation
) # animate map

# Use "here" package to set working directory
here::i_am("healthRiskADAPT.Rproj")

#### LOAD SIMULATION DATA ####

# cleanup
rm(list = ls())

set.seed(42)  # Ensure reproducibility

# get list of EPA stations
if(!file.exists(here("data", "sens_analysis", "grid_data.csv"))) {
    
    # load tract dataset
    df_tract <- read_rds(here("data", "US", "CA", "df_tract.rds")) %>%
        select(!geometry)
    
    # Define the number of simulations
    n_sim <- 25
    
    # Define distributions for house parameters (example)
    # Assume normal distributions for parameters or based on your data
    
    ##### Floor Area #####
    # ggplot(df_tract, aes(x = "", y = house_area)) +
    #     geom_violin(fill = "lightblue", color = "black") +
    #     labs(title = "Violin Plot of House Area",
    #          x = "",
    #          y = "House Area (sq. ft.)") +
    #     theme_minimal()
    
    # ggplot(df_tract, aes(x = "", y = house_area)) +
    #     geom_boxplot(fill = "lightblue", color = "black") +
    #     labs(title = "Violin Plot of House Area",
    #          x = "",
    #          y = "House Area (sq. ft.)") +
    #     theme_minimal()
    
    # ggplot(data.frame(df_tract), aes(x = house_area)) +
    #     geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    #     labs(title = "Monte Carlo Sampled Distribution of House Area",
    #          x = "House Area (sq. ft.)",
    #          y = "Frequency") +
    #     theme_minimal()
    
    ###### MC house area ######
    
    # Generate Monte Carlo samples (e.g., 1,000)
    # mc_house_area_samples <- sample(df_tract$house_area, size = n_sim, replace = TRUE)
    
    # Plot the Monte Carlo Sample Distribution
    # ggplot(data.frame(mc_house_area_samples), aes(x = mc_house_area_samples)) +
    #     geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    #     labs(title = "Monte Carlo Sampled Distribution of House Area",
    #          x = "House Area (sq. ft.)",
    #          y = "Frequency") +
    #     theme_minimal()
    
    # Combine original and MC sample data into a dataframe
    # df_plot <- data.frame(
    #     value = c(df_tract$house_area, mc_house_area_samples),
    #     type = c(rep("Original Data", length(df_tract$house_area)), 
    #              rep("Monte Carlo Samples", length(mc_house_area_samples)))
    # )
    
    # Create the histogram
    # ggplot(df_plot, aes(x = value, y = ..density..,  fill = type)) +
    #     geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
    #     scale_fill_manual(values = c("Original Data" = "blue", "Monte Carlo Samples" = "red")) +
    #     labs(title = "Comparison of House Area: Original vs Monte Carlo Samples",
    #          x = "House Area (sq. ft.)",
    #          #y = "Frequency",
    #          y = "Density", # if y = ..density..,, otherwise "Frequency"
    #          fill = "Dataset") +
    #     theme_minimal()
    
    ###### Sobol House Area ######
    
    # Generate Sobol sequence (between 0 and 1)
    sobol_samples <- sobol(n = n_sim, dim = 1)  
    
    # Scale Sobol samples to match `house_area` range
    house_min <- min(df_tract$house_area, na.rm = TRUE)
    house_max <- max(df_tract$house_area, na.rm = TRUE)
    
    # Convert Sobol samples to closest real `house_area` values using quantiles
    sobol_house_area_samples <- quantile(df_tract$house_area, sobol_samples, type = 8)
    
    # Combine original and Sobol sampled data
    # df_plot <- data.frame(
    #     value = c(df_tract$house_area, sobol_house_area_samples),
    #     type = c(rep("Original Data", length(df_tract$house_area)), 
    #              rep("Sobol Samples", length(sobol_house_area_samples)))
    # )
    
    # Plot Histogram for Comparison
    # ggplot(df_plot, aes(x = value, y = ..density.., fill = type)) +
    #     geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
    #     scale_fill_manual(values = c("Original Data" = "blue", "Sobol Samples" = "green")) +
    #     labs(title = "Comparison of House Area: Original vs Sobol Samples",
    #          x = "House Area (sq. m)",
    #          #y = "Frequency",
    #          y = "Density", # if y = ..density..,, otherwise "Frequency"
    #          fill = "Dataset") +
    #     theme_minimal()
    
    rm(house_max, house_min, sobol_samples)
    
    
    ##### Effective Leakage #####
    
    # Generate Sobol sequence (between 0 and 1)
    sobol_samples <- sobol(n = n_sim, dim = 1)  
    
    # Scale Sobol samples to match `house_height` range
    house_min <- min(df_tract$effect_leakage, na.rm = TRUE)
    house_max <- max(df_tract$effect_leakage, na.rm = TRUE)
    
    # Convert Sobol samples to closest real `house_height` values using quantiles
    sobol_effect_leakage_samples <- quantile(df_tract$effect_leakage[!is.na(df_tract$effect_leakage)], sobol_samples, type = 8)
    
    # Combine original and Sobol sampled data
    #df_plot <- data.frame(
    #    value = c(df_tract$effect_leakage, sobol_effect_leakage_samples),
    #    type = c(rep("Original Data", length(df_tract$effect_leakage)), 
    #             rep("Sobol Samples", length(sobol_effect_leakage_samples)))
    #)
    
    # Plot Histogram for Comparison
    #ggplot(df_plot, aes(x = value, y = ..density.., fill = type)) +
    #    geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
    #    scale_fill_manual(values = c("Original Data" = "blue", "Sobol Samples" = "green")) +
    #    labs(title = "Comparison of House Effective Leakage: Original vs Sobol Samples",
    #         x = "Square m",
    #         #y = "Frequency",
    #         y = "Density", # if y = ..density..,, otherwise "Frequency"
    #         fill = "Dataset") +
    #    theme_minimal()
    
    rm(house_max, house_min, sobol_samples)
    
    
    ##### MERV FILTERS #####
    "
    We can use MERV 8, MERV 10, MERV 13, MERV 16
    "
    
    
    ##### geoid RANDOM SAMPLE n_sim tracts #####
    
    mc_geoid_samples <- sample(df_tract$geoid, size = n_sim)
    
    
    ##### I/O ratio as first guess #####
    "
    We can use 0.25, 0.5, 0.75
    "
    
    
    ##### Combinations #####
    
    grid_data <- expand_grid(c(0.25, 0.5, 0.75), c(8,10,13,16), mc_geoid_samples, sobol_house_area_samples, sobol_effect_leakage_samples)
    
    grid_data <- grid_data %>% 
        mutate(test = row_number())
    
    names(grid_data) <- c("start_IO", "MERV", "geoid", "house_area", "effect_leakage", "test")
    
    grid_data <- grid_data %>%
        select(test, geoid, start_IO, MERV, house_area, effect_leakage)
    
    write_csv(grid_data, file = here("data", "sens_analysis", "grid_data.csv"))
    
} else {
    grid_data <- read_csv(file = here("data", "sens_analysis", "grid_data.csv"))
}


#### simulations from 7 to 13 sept 2020 ####

#### DATA PREP ####

## FILTER ##

# build table of MERV filters
df_filter <- tibble(filter_merv = c(0, 8, 9, 10, 11, 12, 13, 14, 15, 16),
                    filter_eff = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 0.85, 0.9, 0.9, 0.95),
                    filter_pres = c(0, 0.056, 0.056, 0.056, 0.088, 0.088, 0.088, 0.120, 0.120, 0.120)) # pressure loss


## TRACTS ##

# load tract data
df_tract <- read_rds(here("data", "US", "CA", "df_tract.rds"))
df_tract <- df_tract %>% st_drop_geometry()

# heating load for air flow rate##
df_tract <- df_tract %>%
    mutate(heat_load = case_when(climate_ashrae == "2B" ~ 40,
                                 climate_ashrae == "3B" ~ 45,
                                 climate_ashrae == "3C" ~ 45,
                                 climate_ashrae == "4B" ~ 50,
                                 climate_ashrae == "5B" ~ 60,
                                 climate_ashrae == "6B" ~ 65))

# NOTE: (county_new) We approximate the indoor temperature for the 8 counties without Ecobee users
#       using the indoor temperature of homes from a close county.
#       The selection was manual, considering geographical similarities (flatland, mountain range, etc.)
#       when multiple choices were possible.
#       The new variable for county is overwritten
df_tract <- df_tract %>%
    mutate(county = case_when(county == "colusa" | county == "glenn" | county == "tehama" ~ "sutter",
                              county == "del norte" | county == "mendocino" | county == "trinity" ~ "humboldt",
                              county == "mariposa" ~ "tuolumne",
                              county == "modoc" ~ "lassen",
                              TRUE ~ county))

#### Summary stats for the original df_tract ####

# house area
df_tract %>%
    summarise(
        q25 = round(quantile(house_area, 0.25, na.rm = TRUE)),
        q50 = round(quantile(house_area, 0.50, na.rm = TRUE)),
        q75 = round(quantile(house_area, 0.75, na.rm = TRUE))
    )

# house leakage
df_tract %>%
    summarise(
        q25 = round(quantile(effect_leakage, 0.25, na.rm = TRUE), digits = 4),
        q50 = round(quantile(effect_leakage, 0.50, na.rm = TRUE), digits = 4),
        q75 = round(quantile(effect_leakage, 0.75, na.rm = TRUE), digits = 4)
    )  



#### grid_data to create synthetic sample ####
df_tract <- grid_data %>%
    left_join(df_tract %>% select(geoid, 
                                  tract,
                                  county,
                                  #geometry,
                                  lon,
                                  lat,
                                  climate_ca, 
                                  climate_ashrae,
                                  isd_station,
                                  epa_station,
                                  population,
                                  median_income,
                                  house_value,
                                  poverty_perc,
                                  house_total,
                                  house_detached,
                                  house_storey,
                                  house_height,
                                  poverty_status,
                                  heat_load
    ), by = "geoid") %>% glimpse()

df_tract %>% select(house_height) %>% group_by(house_height) %>% summarise(n())
# "
# The house that are 3 mt height are 95%, let's not consider the 5.5 mt
# "
# df_tract <- df_tract %>%
#     filter(house_height == 3)


# calculate norm_leakage, effect_leakage for tracts
df_tract <- df_tract %>%
    mutate(house_volume = house_area * house_height, # redo house volume for tract with missing house area     
           floor_area_ft2 = house_area * 10.764, # convert from m2 to ft2
           pred_heat_load_btu_hr = round((floor_area_ft2*heat_load) / 3000) * 3000, 
           pred_air_flow_rate_cfm_hi = round((pred_heat_load_btu_hr / (60 * ((0.06965 + 0.07935)/2) * ((0.2404 + 0.2401)/2) * 70)) / 50) * 50,
           design_furnace = pred_air_flow_rate_cfm_hi * (0.3048^3)) # convfrt from cfm to cmm

# calculate norm_leakage, effect_leakage for tracts
df_tract <- df_tract %>%
    mutate(house_volume = house_area * house_height, # redo house volume for tract with missing house area     
           floor_area_ft2 = house_area * 10.764, # convert from m2 to ft2
           pred_heat_load_btu_hr = round((floor_area_ft2*heat_load) / 3000) * 3000, 
           pred_air_flow_rate_cfm_hi = round((pred_heat_load_btu_hr / (60 * ((0.06965 + 0.07935)/2) * ((0.2404 + 0.2401)/2) * 70)) / 50) * 50,
           design_furnace = pred_air_flow_rate_cfm_hi * (0.3048^3)) # convfrt from cfm to cmm

#### only keep columns needed for pm2.5 calculation ####
df_tract <- df_tract %>%
    select(test, start_IO, MERV, geoid, county, epa_station, isd_station, house_height, house_volume, house_area, effect_leakage, design_furnace)


#### ISD ####

# load isd dataset
df_isd <- read_rds(here("data", "US", "df_isd.rds"))

## expand to 10 min timesteps and interpolate pm2.5
df_isd <- df_isd %>%
    group_by(isd_station) %>%
    arrange(datetime, .by_group = TRUE) %>%
    complete(datetime = seq(min(datetime), max(datetime), by = "10 min")) %>% #NOTE: find the median on ecobee dataset for 2020 wildfire season
    mutate(wind_speed = na.approx(wind_speed, rule = 2),
           t_out = na.approx(t_out, rule = 2)) %>%
    ungroup()


#### EPA PM2.5 ####

# load epa data
df_epa <- read_rds(here("data", "US", "df_epa.rds")) 

# replace NaNs with interpolated values
df_epa <- df_epa %>%
    mutate(pm_25_out = ifelse(is.nan(pm_25_out), NA, pm_25_out),
           pm_25_out = na.approx(pm_25_out, rule = 2))

## expand to 10 min timesteps and interpolate pm2.5
df_epa <- df_epa %>%
    group_by(epa_station) %>%
    arrange(datetime, .by_group = TRUE) %>%
    complete(datetime = seq(min(datetime), max(datetime), by = "10 min")) %>%
    mutate(pm_25_out = na.approx(pm_25_out, rule = 2)) %>%
    ungroup()


#### ECOBEE ####

# load ecobee data for entire CA
df_ecobee <- read_rds(here("data", "US", "df_ecobee.rds"))

# make averages for date_time and county
df_ecobee <- df_ecobee %>%
    drop_na(county) %>% # remove row with county == NA
    group_by(date_time, county) %>%
    summarise(t_in = mean(t_ctrl, na.rm = TRUE)) %>%
    ungroup()

# expand to 10 mins and interpolate temperatures
df_ecobee <- df_ecobee %>%
    group_by(county) %>%
    arrange(date_time, .by_group = TRUE) %>%
    complete(date_time = seq(min(date_time), max(date_time), by = "10 min")) %>%
    mutate(t_in = na.approx(t_in, rule = 2)) %>%
    ungroup()

# fan runtime for the baseline scenario
df_counties_baseline_run <- read_rds(here("data", "US", "df_counties_baseline_run.rds"))


#### DYNAMIC PM2.5 MODEL ####

# define list of model constants for stack and wind effects
#NOTE# values pulled from https://doi.org/10.1038/jes.2016.49 and ** from supporting information
model_const <- lst(c_g = 9.8, # gravitational acceleration (m2/s)
                   c_r = 0.5, # fraction of the total leackage area (unitless)**
                   c_t0 = 298, # reference temperature (K)
                   c_x = 0.25, # difference between the fraction of the total leakage area that is in the ceiling and the floor (unitless)**
                   c_c = 0.19, # shielding parameter, referring to the wind shielding resulting from obstructions immediately around the house (unitless)**
                   c_factor = 10, # factor (m)
                   c_a = 0.67, # related to the “terrain class”, which characterizes the general neighborhood where the home is situated (unitless)**
                   c_b = 0.25) # related to the “terrain class”, which characterizes the general neighborhood where the home is situated (unitless)**

# make data directory for the sens_analysis simulations
if(!dir.exists(here("data/sens_analysis/sa_simulations/"))) {
    dir.create(here("data/sens_analysis/sa_simulations/"))
}

# check with geoid are already present in the directory to don't redo the simulations
# get list of simulation files
# list_sim_file <- list.files(path = here("data/sens_analysis/sa_simulations/"), pattern = "[0-9].rds", full.names = FALSE) %>%
#     str_replace(".rds", "")
# 
# df_tract_reduced <- df_tract %>%
#     filter(!test %in% list_sim_file)

# focus on specific time
start_time_sens_analysis <- as.POSIXct("2020-09-07 00:00:00")
end_time_sens_analysis   <- as.POSIXct("2020-09-11 00:00:00")


# parallel
plan(multisession, workers = parallel::detectCores() - 1)  # Use available cores

# run loop for each census tract
#for(i in seq_len(nrow(df_tract))){
#for(i in seq_len(5)){
#for(i in seq_len(nrow(df_tract_reduced))){
process_tract <- function(i) { # for parallelization
    
    print(i)
    
    # get tract info
    sim_house <- df_tract %>%
        slice(i)
    
    # print 
    print(sim_house$test)
    
    # use pm2.5 data to build time series
    sim_house <- df_epa %>%
        left_join(sim_house, ., by = "epa_station")
    
    # add indoor temperatures from ecobee dataset
    # NOTE: see note at: (goto "county_new")
    sim_house <- df_ecobee %>%
        left_join(sim_house, ., by = c("county", "datetime" = "date_time")) %>%
        drop_na(t_in)
    
    # add outdoor temperature and windspeed from ISD dataset
    sim_house <- df_isd %>%
        left_join(sim_house, ., by = c("isd_station", "datetime")) %>%
        drop_na(t_out)
    
    # calculate specific infiltration (m/s), aer
    sim_house <- sim_house %>%
        mutate(d_t = as.numeric(difftime(datetime, lag(datetime, n = 1))), # timesteps
               d_t = replace_na(d_t, median(d_t, na.rm = TRUE)),
               stack_effect = ((1 + (model_const$c_r/2) / 3) * (1 - (model_const$c_x^2)/(2 - model_const$c_r)^2)^1.5 * ((model_const$c_g * house_height) / model_const$c_t0)^0.5),
               wind_effect = (model_const$c_c * (1 - model_const$c_r)^0.333) * (model_const$c_a * (house_height / 10)^model_const$c_b),
               s_inf = sqrt(stack_effect^2 * abs(t_out - t_in) + wind_effect^2 * wind_speed^2), # specific infiltration (https://doi.org/10.1038/jes.2016.49, equation nr.7) #FEDE# abs()???
               q_f = effect_leakage * s_inf * 3600, # airflow due to infiltration through small unintentional openings (m3/h) (https://doi.org/10.1038/jes.2016.49, equation nr.3)
               q_nat = 0, # set natural ventilation to 0 (assume closed windows)
               q_tot = sqrt(q_f^2 + q_nat^2),
               aer = q_tot / house_volume, # 1/hr air changes per hour
               v_infil = aer * house_volume / 60 * d_t, # m3, 60 for having aer/min
               house_aspect_ratio = 1.2) # [-] aspect ratio of the house (l/w)
    
    # only keep columns needed for pm2.5 calculations
    sim_house <- sim_house %>% 
        select(test, start_IO, MERV, county, datetime, d_t, house_volume, house_area, v_infil, design_furnace, pm_25_out)
    
    # Focus on specific interval
    sim_house <- sim_house %>%
        filter(datetime >= start_time_sens_analysis & datetime < end_time_sens_analysis)
    
    # add the ecobee county fan runtime
    sim_house <- sim_house %>%
        left_join(., df_counties_baseline_run, by = c("county", "datetime" = "date_time"))
    
    # define model variables 
    sim_house <- sim_house %>%
        mutate(pm25_pene = 0.75, # penetration factor of PM 2.5 particles
               pm25_dep = 0.5, # [- per hour] deposition per hour
               filter_merv = MERV, # merv rating # THIS IS NEW FOR THE SENS ANALYSIS
               n_lr = 0.15, # leakage as a percentage of v_fan at the return side duct
               n_ls = 0.15, # leakage as a percentage of v_fan at the supply side duct
               #pac_cadr = 156, # cfm median CADR for portable air cleaners (source: df_pac.csv)
               pac_cadr = 156 * (0.3048^3), # m3/min
               pac_filt_eff = 0.9997, # FEDE # check this number for EPA filter
               pac_sqft = 242, # median room size sqft raccomandation (source: df_pac.csv)
               pac_mq = round(pac_sqft * 0.092903, digits = 1), # 0.092903 is the conversion factor from sqft to m2
               pac_norm = pac_mq / house_area # normalized area
        ) %>% 
        left_join(., df_filter, by = "filter_merv")
    
    # prepare other values for calculation
    sim_house <- sim_house %>%
        arrange(datetime) %>%
        mutate(
            # global
            v_de = house_volume * pm25_dep/60 * d_t, #m3 
            
            # baseline
            c_in_baseline = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # estimate initial indoor pm2.5 concentration (#FEDE 0.5 instead of 0.35 - be conservative) ### SENSITIVITY ANALYSIS WE CAN USE start_IO
            v_fan_baseline = ifelse(minute(datetime) < 0, yes = 1, no = 0) * design_furnace * d_t, # m3
            v_da_baseline = v_fan_baseline * ((1 - filter_pres) - n_ls),
            v_ra_baseline = (v_infil + v_da_baseline) * as.integer(v_fan_baseline > 0), # m3
            c_bf_baseline = 0, # (v_ra_baseline*c_in_baseline + v_fan_baseline*n_lr*pm_25_out) / (v_ra_baseline + v_fan_baseline*n_lr), 
            c_af_baseline = c_bf_baseline*(1 - filter_eff),
            
            # 10-min hvac 
            c_in_hvac = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration
            v_fan_hvac = ifelse(minute(datetime) < 10, yes = 1, no = 0) * design_furnace * d_t, # m3 ##FEDE remove hard coded 10.. 
            v_da_hvac = v_fan_hvac * ((1 - filter_pres) - n_ls), # m3
            v_ra_hvac = (v_infil + v_da_hvac) * as.integer(v_fan_hvac > 0), # m3, second part evaluating that is non-zero only when HVAC is ON
            c_bf_hvac = case_when( (v_ra_hvac + v_fan_hvac*n_lr) > 0 ~ (v_ra_hvac*c_in_hvac + v_fan_hvac*n_lr*pm_25_out)/(v_ra_hvac + v_fan_hvac*n_lr), # ug/m3
                                   (v_ra_hvac + v_fan_hvac*n_lr) == 0 ~ 0,
                                   TRUE ~ 0
            ),
            c_af_hvac = c_bf_hvac*(1 - filter_eff), # ug/m3
            
            # baseline county hvac NEW
            c_in_base_hvac = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
            v_fan_base_hvac = fan_baseline_logic * design_furnace * d_t, # m3 
            v_da_base_hvac = v_fan_base_hvac * ((1 - filter_pres) - n_ls), # m3
            v_ra_base_hvac = (v_infil + v_da_base_hvac) * as.integer(v_fan_base_hvac > 0), # m3, second part evaluating that is non-zero only when HVAC is ON
            c_bf_base_hvac = case_when( (v_ra_base_hvac + v_fan_base_hvac*n_lr) > 0 ~ (v_ra_base_hvac*c_in_base_hvac + v_fan_base_hvac*n_lr*pm_25_out)/(v_ra_base_hvac + v_fan_base_hvac*n_lr), # ug/m3
                                        (v_ra_base_hvac + v_fan_base_hvac*n_lr) == 0 ~ 0,
                                        TRUE ~ 0
            ),
            c_af_base_hvac = c_bf_base_hvac*(1 - filter_eff), # ug/m3
            
            # control logic hvac with BASELINE and outdoor PM2.5 >= 35
            c_in_base_out_logic = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
            #v_fan_logic = ifelse(pm_25_out > 35, yes = 1, no = 0) * design_furnace * d_t, # m3
            v_fan_base_out_logic = case_when(fan_baseline_logic == 1 ~ 1, # when the HVAC is running for heating/cooling
                                             pm_25_out >= 35 ~ 1,
                                             pm_25_out < 35 ~ 0,
                                             TRUE ~ 0) * design_furnace * d_t, # m3
            v_da_base_out_logic = v_fan_base_out_logic * ((1 - filter_pres) - n_ls), # m3
            v_ra_base_out_logic = (v_infil + v_da_base_out_logic) * as.integer(v_fan_base_out_logic > 0), # m3
            #c_bf_logic = (v_ra_logic*c_in_logic + v_fan_logic*n_lr*pm_25_out)/(v_ra_logic + v_fan_logic*n_lr), # ug/m3 # we have zero division
            c_bf_base_out_logic = case_when((v_ra_base_out_logic + v_fan_base_out_logic*n_lr) > 0 ~ (v_ra_base_out_logic*c_in_base_out_logic + v_fan_base_out_logic*n_lr*pm_25_out)/(v_ra_base_out_logic + v_fan_base_out_logic*n_lr),
                                            (v_ra_base_out_logic + v_fan_base_out_logic*n_lr) == 0 ~ 0,
                                            TRUE ~ 0
            ),
            c_af_base_out_logic = c_bf_base_out_logic*(1 - filter_eff), # ug/m3
            
            # control logic hvac with BASELINE and indoor PM2.5 >= 5
            c_in_base_in_logic = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
            #v_fan_logic = ifelse(pm_25_out > 35, yes = 1, no = 0) * design_furnace * d_t, # m3
            v_fan_base_in_logic = case_when(fan_baseline_logic == 1 ~ 1, # when the HVAC is running for heating/cooling
                                            c_in_base_in_logic >= 5 ~ 1,
                                            c_in_base_in_logic < 5 ~ 0,
                                            TRUE ~ 0) * design_furnace * d_t, # m3
            v_da_base_in_logic = v_fan_base_in_logic * ((1 - filter_pres) - n_ls), # m3
            v_ra_base_in_logic = (v_infil + v_da_base_in_logic) * as.integer(v_fan_base_in_logic > 0), # m3
            #c_bf_logic = (v_ra_logic*c_in_logic + v_fan_logic*n_lr*pm_25_out)/(v_ra_logic + v_fan_logic*n_lr), # ug/m3 # we have zero division
            c_bf_base_in_logic = case_when((v_ra_base_in_logic + v_fan_base_in_logic*n_lr) > 0 ~ (v_ra_base_in_logic*c_in_base_in_logic + v_fan_base_in_logic*n_lr*pm_25_out)/(v_ra_base_in_logic + v_fan_base_in_logic*n_lr),
                                           (v_ra_base_in_logic + v_fan_base_in_logic*n_lr) == 0 ~ 0,
                                           TRUE ~ 0
            ),
            c_af_base_in_logic = c_bf_base_in_logic*(1 - filter_eff), # ug/m3
            
            # control logic hvac NO BASELINE FAN RUNTIME - ONLY AQ MODE
            c_in_logic = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
            #v_fan_logic = ifelse(pm_25_out > 35, yes = 1, no = 0) * design_furnace * d_t, # m3
            v_fan_logic = case_when(pm_25_out >= 35 & c_in_logic >= 5 ~ 1, 
                                    pm_25_out >= 35 & c_in_logic < 5 ~ 0, 
                                    pm_25_out < 35 & c_in_logic >= 5 ~ 1, 
                                    pm_25_out < 35 & c_in_logic < 5 ~ 0,  
                                    TRUE ~ 0
            ) * design_furnace * d_t, # m3
            v_da_logic = v_fan_logic * ((1 - filter_pres) - n_ls), # m3
            v_ra_logic = (v_infil + v_da_logic) * as.integer(v_fan_logic > 0), # m3
            #c_bf_logic = (v_ra_logic*c_in_logic + v_fan_logic*n_lr*pm_25_out)/(v_ra_logic + v_fan_logic*n_lr), # ug/m3 # we have zero division
            c_bf_logic = case_when((v_ra_logic + v_fan_logic*n_lr) > 0 ~ (v_ra_logic*c_in_logic + v_fan_logic*n_lr*pm_25_out)/(v_ra_logic + v_fan_logic*n_lr),
                                   (v_ra_logic + v_fan_logic*n_lr) == 0 ~ 0,
                                   TRUE ~ 0
            ),
            c_af_logic = c_bf_logic*(1 - filter_eff), # ug/m3
            
            # Portable Air Purifier
            v_fan_pac = ifelse(minute(datetime) <= 60 , yes = 1, no = 0) * ( pac_cadr / pac_filt_eff * pac_norm) * d_t, # m3 
            v_da_pac = v_fan_pac, # m3
            v_ra_pac = v_fan_pac, # m3
            c_in_pac = if_else(row_number() == 1, pm_25_out*start_IO, NA_real_), # estimate initial indoor pm2.5 concentration 
        )
    
    # run simulation loop 
    
    {
        tmp_sim_startime <- Sys.time()
        
        for (i in seq_len(nrow(sim_house)-1)) {
            
            # baseline, no HVAC
            sim_house[i+1,"c_bf_baseline"] <- (sim_house$v_ra_baseline[i]*sim_house$c_in_baseline[i] + sim_house$v_fan_baseline[i]*sim_house$n_lr[i]*sim_house$pm_25_out[i]) /
                (sim_house$v_ra_baseline[i] + sim_house$v_fan_baseline[i]*sim_house$n_lr[i])
            sim_house[i+1,"c_bf_baseline"] <- ifelse(is.na(sim_house$c_bf_baseline[i]), 0, sim_house$c_bf_baseline[i])
            sim_house[i+1,"c_af_baseline"] <- sim_house$c_bf_baseline[i]*( 1 - (sim_house$filter_eff[i]) )
            sim_house[i+1,"c_in_baseline"] <- (
                sim_house$c_in_baseline[i]*sim_house$house_volume[i] - sim_house$v_ra_baseline[i]*(sim_house$c_in_baseline[i]) + 
                    sim_house$v_da_baseline[i]*sim_house$c_af_baseline[i] + sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_baseline[i]) - sim_house$v_de[i]*(sim_house$c_in_baseline[i]) 
            ) / sim_house$house_volume[i]
            sim_house[i+1,"c_in_baseline"] <- ifelse(sim_house[i+1,"c_in_baseline"] < 0, 0, sim_house[i+1,"c_in_baseline"])
            
            # HVAC running for fixed mins/hr 
            ### 1- Forecast the indoor concentration at the next step
            sim_house[i+1,"c_in_hvac"] <- (
                sim_house$c_in_hvac[i]*sim_house$house_volume[i] - sim_house$v_ra_hvac[i]*(sim_house$c_in_hvac[i]) + 
                    sim_house$v_da_hvac[i]*sim_house$c_af_hvac[i] + sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_hvac[i]) - sim_house$v_de[i]*(sim_house$c_in_hvac[i]) 
            ) / sim_house$house_volume[i] # ug/m3
            sim_house[i+1,"c_in_hvac"] <- ifelse(sim_house[i+1,"c_in_hvac"] < 0, 0, sim_house[i+1,"c_in_hvac"]) # NOTE.. is the deposition making the value of concentration going below zero?
            ### 2- Calculate concentration on ducts *Fan runtime is known and fixed as well as V_da and V_ra
            sim_house[i+1,"c_bf_hvac"] <- case_when( (sim_house$v_ra_hvac[i+1] + sim_house$v_fan_hvac[i+1]*sim_house$n_lr[i+1]) > 0 ~ (sim_house$v_ra_hvac[i+1]*sim_house$c_in_hvac[i+1] + sim_house$v_fan_hvac[i+1]*sim_house$n_lr[i+1]*sim_house$pm_25_out[i+1]) / 
                                                         (sim_house$v_ra_hvac[i+1] + sim_house$v_fan_hvac[i+1]*sim_house$n_lr[i+1]),
                                                     (sim_house$v_ra_hvac[i+1] + sim_house$v_fan_hvac[i+1]*sim_house$n_lr[i+1]) == 0 ~ 0,
                                                     TRUE ~ 0
            )
            sim_house[i+1,"c_af_hvac"] <- sim_house$c_bf_hvac[i+1]*( 1 - (sim_house$filter_eff[i+1]) )
            
            # HVAC running using the average run time of the County
            ### 1- Forecast the indoor concentration at the next step
            sim_house[i+1,"c_in_base_hvac"] <- (
                sim_house$c_in_base_hvac[i]*sim_house$house_volume[i] - sim_house$v_ra_base_hvac[i]*(sim_house$c_in_base_hvac[i]) + 
                    sim_house$v_da_base_hvac[i]*sim_house$c_af_base_hvac[i] + sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_base_hvac[i]) - sim_house$v_de[i]*(sim_house$c_in_base_hvac[i]) 
            ) / sim_house$house_volume[i] # ug/m3
            sim_house[i+1,"c_in_base_hvac"] <- ifelse(sim_house[i+1,"c_in_base_hvac"] < 0, 0, sim_house[i+1,"c_in_base_hvac"]) # Correct if the deposition is making the value of concentration going below zero?
            ### 2- Calculate concentration on ducts *Fan runtime is known and fixed as well as V_da and V_ra
            sim_house[i+1,"c_bf_base_hvac"] <- case_when( (sim_house$v_ra_base_hvac[i+1] + sim_house$v_fan_base_hvac[i+1]*sim_house$n_lr[i+1]) > 0 ~ (sim_house$v_ra_base_hvac[i+1]*sim_house$c_in_base_hvac[i+1] + sim_house$v_fan_base_hvac[i+1]*sim_house$n_lr[i+1]*sim_house$pm_25_out[i+1]) / 
                                                              (sim_house$v_ra_base_hvac[i+1] + sim_house$v_fan_base_hvac[i+1]*sim_house$n_lr[i+1]),
                                                          (sim_house$v_ra_base_hvac[i+1] + sim_house$v_fan_base_hvac[i+1]*sim_house$n_lr[i+1]) == 0 ~ 0,
                                                          TRUE ~ 0
            )
            sim_house[i+1,"c_af_base_hvac"] <- sim_house$c_bf_base_hvac[i+1]*( 1 - (sim_house$filter_eff[i+1]) )
            
            
            # HVAC control logic over baseline and when outdoor PM2.5 >= 35 ug/m3
            ### 1- Forecast the indoor concentration at the next step
            sim_house[i+1,"c_in_base_out_logic"] <- (
                sim_house$c_in_base_out_logic[i]*sim_house$house_volume[i] - sim_house$v_ra_base_out_logic[i]*(sim_house$c_in_base_out_logic[i]) + 
                    sim_house$v_da_base_out_logic[i]*sim_house$c_af_base_out_logic[i] + sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_base_out_logic[i]) - sim_house$v_de[i]*(sim_house$c_in_base_out_logic[i]) 
            ) / sim_house$house_volume[i] # ug/m3
            sim_house[i+1,"c_in_base_out_logic"] <- ifelse(sim_house[i+1,"c_in_base_out_logic"] < 0, 0, sim_house[i+1,"c_in_base_out_logic"])
            ### 2- Decide the logic to run based on indoor thermal comfort outdoor PM2.5 concentration
            sim_house[i+1, "v_fan_base_out_logic"] <- case_when(sim_house$v_fan_base_out_logic[i+1] > 0 ~ 1, # here there is an optimization possible as we are overwriting an already present calculation
                                                                sim_house$pm_25_out[i+1] >= 35 ~ 1,
                                                                sim_house$pm_25_out[i+1] <  35 ~ 0,
                                                                TRUE ~ 0) * sim_house$design_furnace[i+1] * sim_house$d_t[i+1] # m3
            ### 3- Calculate flows and concentrations
            sim_house[i+1, "v_da_base_out_logic"] <- sim_house$v_fan_base_out_logic[i+1] * ((1 - sim_house$filter_pres[i+1]) - sim_house$n_ls[i+1])
            sim_house[i+1, "v_ra_base_out_logic"] <- (sim_house$v_infil[i+1] + sim_house$v_da_base_out_logic[i+1]) * as.integer(sim_house$v_fan_base_out_logic[i+1] > 0)
            sim_house[i+1,"c_bf_base_out_logic"] <- case_when( (sim_house$v_ra_base_out_logic[i+1] + sim_house$v_fan_base_out_logic[i+1]*sim_house$n_lr[i+1]) > 0 ~ (sim_house$v_ra_base_out_logic[i+1]*sim_house$c_in_base_out_logic[i+1] + sim_house$v_fan_base_out_logic[i+1]*sim_house$n_lr[i+1]*sim_house$pm_25_out[i+1]) / 
                                                                   (sim_house$v_ra_base_out_logic[i+1] + sim_house$v_fan_base_out_logic[i+1]*sim_house$n_lr[i+1]),
                                                               (sim_house$v_ra_base_out_logic[i+1] + sim_house$v_fan_base_out_logic[i+1]*sim_house$n_lr[i+1]) == 0 ~ 0,
                                                               TRUE ~ 0
            )
            sim_house[i+1,"c_af_base_out_logic"] <- sim_house$c_bf_base_out_logic[i+1]*( 1 - (sim_house$filter_eff[i+1]) )
            
            # HVAC control logic over baseline and when indoor PM2.5 >= 5 ug/m3
            ### 1- Forecast the indoor concentration at the next step
            sim_house[i+1,"c_in_base_in_logic"] <- (
                sim_house$c_in_base_in_logic[i]*sim_house$house_volume[i] - sim_house$v_ra_base_in_logic[i]*(sim_house$c_in_base_in_logic[i]) + 
                    sim_house$v_da_base_in_logic[i]*sim_house$c_af_base_in_logic[i] + sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_base_in_logic[i]) - sim_house$v_de[i]*(sim_house$c_in_base_in_logic[i]) 
            ) / sim_house$house_volume[i] # ug/m3
            sim_house[i+1,"c_in_base_in_logic"] <- ifelse(sim_house[i+1,"c_in_base_in_logic"] < 0, 0, sim_house[i+1,"c_in_base_in_logic"])
            ### 2- Decide the logic to run based on indoor thermal comfort indoor PM2.5 concentration
            sim_house[i+1, "v_fan_base_in_logic"] <- case_when(sim_house$v_fan_base_in_logic[i+1] > 0 ~ 1, # here there is an optimization possible as we are overwriting an already present calculation
                                                               sim_house$c_in_base_in_logic[i+1] >= 5 ~ 1,
                                                               sim_house$c_in_base_in_logic[i+1] <  5 ~ 0,
                                                               TRUE ~ 0) * sim_house$design_furnace[i+1] * sim_house$d_t[i+1] # m3
            ### 3- Calculate flows and concentrations
            sim_house[i+1, "v_da_base_in_logic"] <- sim_house$v_fan_base_in_logic[i+1] * ((1 - sim_house$filter_pres[i+1]) - sim_house$n_ls[i+1])
            sim_house[i+1, "v_ra_base_in_logic"] <- (sim_house$v_infil[i+1] + sim_house$v_da_base_in_logic[i+1]) * as.integer(sim_house$v_fan_base_in_logic[i+1] > 0)
            sim_house[i+1,"c_bf_base_in_logic"] <- case_when( (sim_house$v_ra_base_in_logic[i+1] + sim_house$v_fan_base_in_logic[i+1]*sim_house$n_lr[i+1]) > 0 ~ (sim_house$v_ra_base_in_logic[i+1]*sim_house$c_in_base_in_logic[i+1] + sim_house$v_fan_base_in_logic[i+1]*sim_house$n_lr[i+1]*sim_house$pm_25_out[i+1]) / 
                                                                  (sim_house$v_ra_base_in_logic[i+1] + sim_house$v_fan_base_in_logic[i+1]*sim_house$n_lr[i+1]),
                                                              (sim_house$v_ra_base_in_logic[i+1] + sim_house$v_fan_base_in_logic[i+1]*sim_house$n_lr[i+1]) == 0 ~ 0,
                                                              TRUE ~ 0
            )
            sim_house[i+1,"c_af_base_in_logic"] <- sim_house$c_bf_base_in_logic[i+1]*( 1 - (sim_house$filter_eff[i+1]) )
            
            
            # HVAC control logic - AQ mode, no HVAC baseline
            ### 1- Forecast the indoor concentration at the next step
            sim_house[i+1,"c_in_logic"] <- (
                sim_house$c_in_logic[i]*sim_house$house_volume[i] - sim_house$v_ra_logic[i]*(sim_house$c_in_logic[i]) + 
                    sim_house$v_da_logic[i]*sim_house$c_af_logic[i] + sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_logic[i]) - sim_house$v_de[i]*(sim_house$c_in_logic[i]) 
            ) / sim_house$house_volume[i] # ug/m3
            sim_house[i+1,"c_in_logic"] <- ifelse(sim_house[i+1,"c_in_logic"] < 0, 0, sim_house[i+1,"c_in_logic"])
            ### 2- Decide the logic to run based on outdoor and indoor concentration
            sim_house[i+1, "v_fan_logic"] <- case_when(sim_house$pm_25_out[i+1] >= 35 & sim_house$c_in_logic[i+1] >= 5 ~ 1, 
                                                       sim_house$pm_25_out[i+1] >= 35 & sim_house$c_in_logic[i+1] < 5 ~ 0,  
                                                       sim_house$pm_25_out[i+1] <  35 & sim_house$c_in_logic[i+1] >= 5 ~ 1,  
                                                       sim_house$pm_25_out[i+1] <  35 & sim_house$c_in_logic[i+1] < 5 ~ 0,   
                                                       TRUE ~ 0) * sim_house$design_furnace[i+1] * sim_house$d_t[i+1] # m3
            ### 3- Calculate flows and concentrations
            sim_house[i+1, "v_da_logic"] <- sim_house$v_fan_logic[i+1] * ((1 - sim_house$filter_pres[i+1]) - sim_house$n_ls[i+1])
            sim_house[i+1, "v_ra_logic"] <- (sim_house$v_infil[i+1] + sim_house$v_da_logic[i+1]) * as.integer(sim_house$v_fan_logic[i+1] > 0)
            sim_house[i+1,"c_bf_logic"] <- case_when( (sim_house$v_ra_logic[i+1] + sim_house$v_fan_logic[i+1]*sim_house$n_lr[i+1]) > 0 ~ (sim_house$v_ra_logic[i+1]*sim_house$c_in_logic[i+1] + sim_house$v_fan_logic[i+1]*sim_house$n_lr[i+1]*sim_house$pm_25_out[i+1]) / 
                                                          (sim_house$v_ra_logic[i+1] + sim_house$v_fan_logic[i+1]*sim_house$n_lr[i+1]),
                                                      (sim_house$v_ra_logic[i+1] + sim_house$v_fan_logic[i+1]*sim_house$n_lr[i+1]) == 0 ~ 0,
                                                      TRUE ~ 0
            )
            sim_house[i+1,"c_af_logic"] <- sim_house$c_bf_logic[i+1]*( 1 - (sim_house$filter_eff[i+1]) )
            
            # PAC running all the time
            sim_house[i+1,"c_in_pac"] <- (
                sim_house$c_in_pac[i]*sim_house$house_volume[i] - # ug, tot numb particle house
                    sim_house$v_fan_pac[i]*(sim_house$c_in_pac[i]) * (sim_house$pac_filt_eff[i]) + # ug, total number particles removed by the filter
                    sim_house$v_infil[i]*sim_house$pm_25_out[i]*sim_house$pm25_pene[i] - # ug, particle infiltration, # NOTE 75% particle go trough the walls
                    sim_house$v_infil[i]*(sim_house$c_in_pac[i]) - # ug, particle outfiltration
                    sim_house$v_de[i]*(sim_house$c_in_pac[i]) # ug, particle deposition
            ) / sim_house$house_volume[i]
            sim_house[i+1,"c_in_pac"] <- ifelse(sim_house[i+1,"c_in_pac"] < 0, 0, sim_house[i+1,"c_in_pac"])
        }
        
        tmp_sim_endtime <- Sys.time()
        
        print(tmp_sim_endtime - tmp_sim_startime)
        }
    
    # [DEBUG] plot the results
    # sim_house %>%
    #   filter(datetime > "2020-09-05 00:00:00" & datetime < "2020-09-07 00:00:00" ) %>%
    # ggplot(data = ., aes(x = datetime)) + 
    #   geom_hline(yintercept = 5, lty=2, alpha = 0.6) +
    #   geom_hline(yintercept = 35, lty=2, alpha = 0.6) +
    #   geom_line(aes(y = pm_25_out), alpha = 0.9) + 
    #   geom_line(aes(y=c_in_baseline), color="blue", alpha=0.5) + 
    #   geom_line(aes(y=c_in_hvac), color="firebrick", alpha=0.5) + 
    #   geom_line(aes(y=c_in_logic), color="forestgreen") + 
    #   geom_line(aes(y=c_in_pac), color="purple", alpha=0.5) +
    #   ggtitle(paste0("Outdoor/indoor scenario PM2.5 concentration\nGeoid: ", sim_house$geoid[1], ", ", round(sim_house$house_volume[1], digits = 1), " m3, ", 
    #                  round(sim_house$house_area[1], digits = 1), " m2, ", 
    #                  sim_house$house_volume[1] / sim_house$house_area[1], " h\n", 
    #                  sim_house$datetime[1], " - ",
    #                  sim_house$datetime[nrow(sim_house)]),
    #           subtitle = "Black: outdoor PM2.5, blue: no HVAC, red: HVAC 12 min/hr, green: HVAC control logic, purple: PAC 60 min/hr
    #           Indoor PM2.5 threshold is 5 ug/m3
    #           Outdoor PM2.5 threshold is 35 ug/m3") +
    #   theme(axis.title.x = element_blank(),
    #         axis.title.y = element_blank())
    
    # select columns to keep -> rename to be consistent with the new names 2023-05-27
    sim_house <- sim_house %>%
        mutate(v_fan_hvac = case_when(v_fan_hvac > 0 ~ 1, TRUE ~ 0),
               v_fan_base_hvac = case_when(v_fan_base_hvac > 0 ~ 1, TRUE ~ 0),
               v_fan_base_out_logic = case_when(v_fan_base_out_logic > 0 ~ 1, TRUE ~ 0),
               v_fan_base_in_logic = case_when(v_fan_base_in_logic > 0 ~ 1, TRUE ~ 0),
               v_fan_logic = case_when(v_fan_logic > 0 ~ 1, TRUE ~ 0),
               v_fan_pac = case_when(v_fan_pac > 0 ~ 1, TRUE ~ 0),
        ) %>%
        select(test, datetime, pm_25_out,
               # save infiltration data
               #fan_baseline = v_fan_baseline,
               pm_25_infiltration = c_in_baseline, 
               # save 10-min fixed HVAC runtime 
               fan_10min_hvac = v_fan_hvac,
               pm_25_10min_hvac = c_in_hvac,
               # save baseline HVAC runtime
               fan_baseline_hvac = v_fan_base_hvac,
               pm_25_baseline_hvac = c_in_base_hvac,
               # save baseline and outdoor PM2.5 logic
               fan_baseline_out = v_fan_base_out_logic,
               pm_25_baseline_out = c_in_base_out_logic,
               # save baseline and indoor PM2.5 logic
               fan_baseline_in = v_fan_base_in_logic,
               pm_25_baseline_in = c_in_base_in_logic,
               # save AQ mode and PAC
               fan_AQ_logic = v_fan_logic,
               pm_25_AQ_logic = c_in_logic,
               fan_pac = v_fan_pac,
               pm_25_pac = c_in_pac
        )
    
    # save result
    write_rds(sim_house, file = here("data", "sens_analysis", "sa_simulations", paste0(sim_house[1,1],'.rds')))
    
    # summarise exposure and write to file # FEDE check runtime
    sim_house %>%
        summarise(test = first(test),
                  sample = n(),
                  # save infiltration PM2.5
                  pm_25_infiltration = sum(pm_25_infiltration, na.rm = TRUE),
                  # save total runtime and indoor PM2.5 for the 10-min logic
                  fan_10min_hvac = sum(fan_10min_hvac, na.rm = TRUE) * 10, # to have the total min
                  pm_25_10min_hvac = sum(pm_25_10min_hvac, na.rm = TRUE),
                  # save total runtime and indoor PM2.5 for the heating and air conditioning baseline
                  fan_baseline_hvac = sum(fan_baseline_hvac, na.rm = TRUE) * 10, 
                  pm_25_baseline_hvac = sum(pm_25_baseline_hvac, na.rm = TRUE),
                  # save total runtime and indoor PM2.5 for the baseline scenario and the outdoor PM2.5 control logic
                  fan_baseline_out = sum(fan_baseline_out, na.rm = TRUE) * 10,
                  pm_25_baseline_out = sum(pm_25_baseline_out, na.rm = TRUE),
                  # save total runtime and indoor PM2.5 for the baseline scenario and the indoor PM2.5 control logic
                  fan_baseline_in = sum(fan_baseline_in, na.rm = TRUE) * 10,
                  pm_25_baseline_in = sum(pm_25_baseline_in, na.rm = TRUE),
                  # save total runtime and indoor PM2.5 for the AQ scenario considering only the indoor PM2.5 logic
                  fan_AQ_logic = sum(fan_AQ_logic, na.rm = TRUE) * 10,
                  pm_25_AQ_logic = sum(pm_25_AQ_logic, na.rm = TRUE),
                  # save total runtime and indoor PM2.5 for the PAC
                  fan_pac = sum(fan_pac, na.rm = TRUE) * 10,
                  pm_25_pac = sum(pm_25_pac, na.rm = TRUE),
                  # save the total outdoor PM2.5
                  pm_25_out = sum(pm_25_out, na.rm = TRUE)
        ) %>%
        write_delim(., here("data", "sens_analysis", "sa_simulations", "_total_exposure.csv"), delim = ",", append = TRUE)
    
    return(sim_house)
}

# try running first 5 tests
#results <- future_map(seq_len(5), process_tract, .progress = TRUE)

# try running random 100 test
#random_numbers <- sample(1:180000, 100, replace = FALSE)


# START SIMULATIONS *** CAREFUL ***
results <- future_map(seq_len(nrow(grid_data)), process_tract, .progress = TRUE)

stop("Finishing simulations. Exit source.")


#### MANUAL ####

#### IMPORT DATA ####

# get list of simulation files
list_file <- list.files(path = here("data", "sens_analysis", "sa_simulations"), pattern = "[0-9].rds", full.names = TRUE)

# function to import and summarise data files
import_sim <- function(sim_file) {
    
    # read rds and average over time interval (day, hour etc.)
    df_sim_load <- read_rds(file = sim_file) %>%
        group_by(datetime = floor_date(datetime, unit = "hour")) %>%
        summarise(test = first(test),
                  # infiltration PM2.5
                  pm_25_infiltration = sum(pm_25_infiltration, na.rm = TRUE),
                  # total runtime and indoor PM2.5 for the 10-min logic
                  fan_10min_hvac = sum(fan_10min_hvac, na.rm = TRUE) * 10, # to have the total min
                  pm_25_10min_hvac = sum(pm_25_10min_hvac, na.rm = TRUE),
                  # total runtime and indoor PM2.5 for the heating and air conditioning baseline
                  fan_baseline_hvac = sum(fan_baseline_hvac, na.rm = TRUE) * 10, 
                  pm_25_baseline_hvac = sum(pm_25_baseline_hvac, na.rm = TRUE),
                  # total runtime and indoor PM2.5 for the baseline scenario and the outdoor PM2.5 control logic
                  fan_baseline_out = sum(fan_baseline_out, na.rm = TRUE) * 10,
                  pm_25_baseline_out = sum(pm_25_baseline_out, na.rm = TRUE),
                  # total runtime and indoor PM2.5 for the baseline scenario and the indoor PM2.5 control logic
                  fan_baseline_in = sum(fan_baseline_in, na.rm = TRUE) * 10,
                  pm_25_baseline_in = sum(pm_25_baseline_in, na.rm = TRUE),
                  # total runtime and indoor PM2.5 for the AQ scenario considering only the indoor PM2.5 logic
                  fan_AQ_logic = sum(fan_AQ_logic, na.rm = TRUE) * 10,
                  pm_25_AQ_logic = sum(pm_25_AQ_logic, na.rm = TRUE),
                  # total runtime and indoor PM2.5 for the PAC
                  fan_pac = sum(fan_pac, na.rm = TRUE) * 10,
                  pm_25_pac = sum(pm_25_pac, na.rm = TRUE),
                  # total outdoor PM2.5
                  pm_25_out = sum(pm_25_out, na.rm = TRUE)
        )
    
    # keep dataframe
    return(df_sim_load)
    
}

# run on list of files
list_sim <- pblapply(list_file, import_sim)

# combine data files
df_sim <- bind_rows(list_sim)

# save dataset 
write_rds(df_sim, here("data", "sens_analysis", "sa_simulations", paste0("df_hourly_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")), compress = "gz")

# clean up
rm(list_sim, list_file, import_sim)



#### PLOTS FOR SENSITIVITY ANALYSIS ####

grid_data <- read_csv(file = here("data", "sens_analysis", "grid_data.csv"))

df_sim_plt <- df_sim %>%
    left_join(y = grid_data, by = "test")

##### SUMMARY STATS #####

small_house <- 90
medium_house <- 120
large_house <- 152

airtight_house <- 0.0273
med.leaky_house <- 0.0359 
leaky_house <- 0.0454

df_sim_plt %>% nrow()

# house area
df_sim_plt %>%
    group_by(test) %>% 
    summarise(start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage)) %>%
    select(house_area) %>%
    mutate(house_size = case_when(house_area <= small_house ~ "small",
                                  house_area > small_house & house_area <= large_house ~ "medium",
                                  house_area > large_house ~ "large")) %>%
    select(house_size) %>% group_by(house_size) %>% summarise(n())

# leakage
df_sim_plt %>%
    group_by(test) %>% 
    summarise(start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage)) %>%
    select(effect_leakage) %>%
    mutate(house_leak = case_when(effect_leakage <= airtight_house ~ "airtight",
                                  effect_leakage > airtight_house & effect_leakage <= leaky_house ~ "med. leaky",
                                  effect_leakage > leaky_house ~ "very leaky"
    )) %>% select(house_leak) %>% group_by(house_leak) %>% summarise(n())

# I/O
df_sim_plt %>%
    group_by(test) %>% 
    summarise(start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage)) %>%
    select(start_IO) %>% 
    group_by(start_IO) %>% 
    summarise(n())



##### PLOTS #####

###### Upgrading Filter and using wildfire mode ######

# add variables "house_size" and "house_leak" to dataset
df_sim_round <- df_sim_plt %>% 
    mutate(house_size = case_when(house_area <= small_house ~ "small",
                                  house_area > small_house & house_area <= large_house ~ "medium",
                                  house_area > large_house ~ "large")) %>%
    mutate(house_leak = case_when(effect_leakage <= airtight_house ~ "airtight",
                                  effect_leakage > airtight_house & effect_leakage <= leaky_house ~ "med.leaky",
                                  effect_leakage > leaky_house ~ "very.leaky"
    ))

# calculate the reduction of indoor exposure by upgrading filter and using wildfire mode - against outdoor exposure

df_sim_round %>%
    group_by(MERV) %>%
    summarise(exp_outdoor = sum(pm_25_out),
              exp_sheltering = sum(pm_25_infiltration),
              exp_wildfire_mode = sum(pm_25_baseline_out),
              exp_diff_shelter = 100 - exp_sheltering/exp_outdoor * 100,
              exp_diff_wildfire_mode = 100 - exp_wildfire_mode/exp_outdoor * 100,
              exp_diff = (exp_diff_wildfire_mode - exp_diff_shelter)
    )

# calculate the reduction of indoor exposure by house size

df_sim_round %>%
    #filter(MERV > 8) %>%
    group_by(house_size) %>%
    summarise(exp_outdoor = sum(pm_25_out),
              exp_sheltering = sum(pm_25_infiltration),
              exp_diff_shelter = 100 - exp_sheltering/exp_outdoor * 100
    )

# calculate the reduction of indoor exposure by house leakage

df_sim_round %>%
    filter(MERV > 8) %>%
    group_by(house_leak) %>%
    summarise(exp_outdoor = sum(pm_25_out),
              exp_sheltering = sum(pm_25_infiltration),
              exp_diff_shelter = 100 - exp_sheltering/exp_outdoor * 100
    )





###### IS I/O IMPORTANT ? ######

# prepare dataset to be used for plots
prep_plot_io <- df_sim_plt %>%
    group_by(test) %>% 
    mutate(house_size = case_when(house_area <= small_house ~ "small",
                                  house_area > small_house & house_area <= large_house ~ "medium",
                                  house_area > large_house ~ "large")) %>%
    mutate(house_leak = case_when(effect_leakage <= airtight_house ~ "airtight",
                                  effect_leakage > airtight_house & effect_leakage <= leaky_house ~ "med.leaky",
                                  effect_leakage > leaky_house ~ "very.leaky"
    )) %>%
    summarise(pm_25_infiltration = sum(pm_25_infiltration),
              fan_10min_hvac = sum(fan_10min_hvac),
              pm_25_10min_hvac = sum(pm_25_10min_hvac),
              fan_baseline_hvac = sum(fan_baseline_hvac),
              pm_25_baseline_hvac = sum(pm_25_baseline_hvac),
              fan_baseline_out = sum(fan_baseline_out),
              pm_25_baseline_out = sum(pm_25_baseline_out),
              fan_baseline_in = sum(fan_baseline_in),
              pm_25_baseline_in = sum(pm_25_baseline_in),
              fan_AQ_logic = sum(fan_AQ_logic),
              pm_25_AQ_logic = sum(pm_25_AQ_logic),
              fan_pac = sum(fan_pac),
              pm_25_pac = sum(pm_25_pac),
              pm_25_out = sum(pm_25_out),
              geoid = first(geoid),
              start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage),
              house_size = first(house_size),
              house_leak = first(house_leak)
    ) %>% 
    mutate(exp_red_shelter = 100 - pm_25_infiltration / pm_25_out * 100,
           exp_red_10_min = 100 - pm_25_10min_hvac / pm_25_out * 100,
           exp_red_baseline_hvac = 100 - pm_25_baseline_hvac / pm_25_out * 100,
           exp_red_baseline_out = 100 - pm_25_baseline_out / pm_25_out * 100,
           exp_red_baseline_in = 100 - pm_25_baseline_in / pm_25_out * 100,
           exp_red_baseline_AQ = 100 - pm_25_AQ_logic / pm_25_out * 100,
           exp_red_baseline_pac = 100 - pm_25_pac / pm_25_out * 100,
           # compared to infiltration
           exp_red_hvac_shelter = 100 - pm_25_baseline_hvac / pm_25_infiltration * 100,
           exp_red_wildfire_shelter = 100 - pm_25_baseline_out / pm_25_infiltration * 100,
           exp_red_pac_shelter = 100 - pm_25_pac / pm_25_infiltration * 100
    ) %>% glimpse()


# THIS IS A NICE COMPREHENSIVE PLOT
prep_plot_io %>% 
    ggplot(aes(x = factor(start_IO), fill = factor(house_leak))) +
    # Violin plot grouped by house_leak
    geom_violin(aes(y = exp_red_baseline_out), alpha = 0.5) +
    # Facet by house_size
    facet_wrap(~ factor(house_size)) +
    facet_wrap(~ factor(MERV)) +
    labs(title = "Wildfire Mode and Buildings Characteristics",
         subtitle = "MERV filter (8-16), House Size and Leakage, I/O PM2.5 ratio",
         y = "Exposure reduction compared to being outdoors (%)",
         x = "Initial Indoor/Outdoor PM2.5 concentration"
    ) +
    # Theme for better readability
    theme_minimal()

ggsave(filename = here("products/gfx/sens_analysis_IO_no_points.png"), width = 3600, height = 1800, dpi = 300, units = "px")

# adding points
prep_plot_io %>% 
    ggplot(aes(x = interaction(start_IO, house_leak), fill = factor(house_leak))) +  # Ensure x-axis grouping
    # Violin plot grouped by house_leak and start_IO
    geom_violin(aes(y = exp_red_baseline_out), alpha = 0.5) +
    # Facet by house_size (rows) and MERV (columns)
    facet_grid(rows = vars(house_size), cols = vars(MERV)) +
    # Add points, ensuring they align with the correct violin
    geom_point(aes(y = exp_red_baseline_out,
                   color = factor(house_leak),
                   group = interaction(start_IO, house_leak, house_size, MERV)), 
               alpha = 0.01, 
               position = position_jitterdodge(jitter.width = 1, jitter.height = 0.2, dodge.width = 0.8), 
               size = 0.5) +
    # Manually set colors for house_leak (violin plots)
    scale_fill_manual(values = c(
        "airtight" = "#2E86C1",  
        "med.leaky" = "#F39C12",
        "very.leaky" = "#E74C3C"
    )) +
    scale_color_manual(values = c(
        "airtight" = "#2E86C1",  
        "med.leaky" = "#F39C12",
        "very.leaky" = "#E74C3C"
    )) +
    labs(title = "Wildfire Mode and Buildings Characteristics",
         subtitle = "MERV filter (8-16), House Size and Leakage, I/O PM2.5 ratio",
         y = "Exposure reduction compared to being outdoors (%)",
         x = "Interaction between I/O PM2.5 and House Leakage"
    ) +
    # Improve legend
    guides(
        #color = guide_legend(override.aes = list(size = 1, alpha = 1)),  # Ensure legend points are small & opaque
        color = "none",  # Ensure legend points are small & opaque
        fill = guide_legend(title = "House Leakiness")  # Custom title for fill legend
    ) +
    # Theme for better readability
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x labels
    )

ggsave(filename = here("products/gfx/sens_analysis_IO.png"), width = 3600, height = 1800, dpi = 300, units = "px")



###### t.test for start IO ######

tmp_data <- prep_plot_io %>%
    #filter(house_size == "medium", house_leak == "med.leaky", MERV == 16)
    filter(house_size == "medium", house_leak == "med.leaky")

# Split data into two groups: start_IO = 0.25 vs. 0.75
group_025 <- tmp_data %>% filter(start_IO == 0.25) %>% pull(exp_red_baseline_out)
group_075 <- tmp_data %>% filter(start_IO == 0.75) %>% pull(exp_red_baseline_out)

# Perform unpaired t-test (Welch's t-test by default)
t_test_result <- t.test(group_025, group_075, var.equal = FALSE)

# Print results
print(t_test_result)

# Compute Cohen's d
cohen_d_value <- effsize::cohen.d(group_025, group_075)
print(cohen_d_value)

prep_plot_io %>%
    #filter(house_size == "medium", house_leak == "med.leaky", MERV == 16) %>%
    ggplot(aes(x = factor(start_IO), y = exp_red_baseline_out, fill = factor(start_IO))) +
    geom_jitter(aes(color = factor(start_IO)), position = position_jitter(width = 0.3, height = 0), size = 0.1, alpha = 0.05) +  # Jitter points
    geom_violin(alpha = 0.5) +  # Violin plot for distribution
    geom_boxplot(width = 0.2, alpha = 0.7) +  # Boxplot inside violin
    labs(title = "Comparison of Wildfire Mode Exposure Reduction by different I/O",
         x = "Start IO", y = "Exposure Reduction Baseline (%)") +
    theme_minimal()

ggsave(filename = here("products/gfx/sens_analysis_IO_ttest.png"), width = 3600, height = 1800, dpi = 300, units = "px")


###### Figure A.8 New ######

# Select relevant columns and reshape data for plotting
plot_data <- prep_plot_io %>%
    select(starts_with("exp_red"), house_size, MERV) %>%
    pivot_longer(cols = starts_with("exp_red"), names_to = "intervention", values_to = "exp_red") %>%
    filter(intervention %in% c("exp_red_shelter", "exp_red_baseline_hvac", "exp_red_baseline_out", "exp_red_baseline_pac")) %>%
    mutate(
        intervention = factor(intervention, 
                              levels = c("exp_red_shelter", "exp_red_baseline_hvac",
                                         "exp_red_baseline_out", 
                                         "exp_red_baseline_pac"),
                              labels = c("Sheltering indoors", "Air conditioning", "Wildfire mode",
                                         "Portable air cleaner")),
        house_size = factor(house_size, levels = c("small", "medium", "large")),
        MERV = factor(MERV, levels = c("8", "10", "13", "16"))
    )

# Create the boxplot
plot_data %>%
    filter(as.numeric(as.character(MERV)) %in% c(10, 13, 16)) %>%
    ggplot(., aes(x = intervention, y = exp_red, fill = interaction(MERV, house_size) )) +
    #geom_jitter(width = 0.2, size = 0.1, alpha = 0.05) +  # Add jitter for data distribution
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplots without outliers for clarity
    labs(
        title = "Indoor PM2.5 exposure reduction (%) relative to being outdoors",
        subtitle = "Scenarious grouped per house floor area and by type of filter",
        x = "Intervention Type",
        y = "Exposure Reduction (%)",
        fill = "Central air system filter type and house size"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

ggsave(filename = here("products/gfx/figureA8.png"), width = 3600, height = 1800, dpi = 300, units = "px")

# quantify uncertainties of exp_red
plot_data %>%
    filter(as.numeric(as.character(MERV)) %in% c(10, 13, 16)) %>%
    group_by(intervention) %>%
    summarise(mean_red = mean(exp_red),
              sd_red = sd(exp_red)
    )

#quantify uncertainties of upgrading MERV
plot_data %>%
    filter(as.numeric(as.character(MERV)) == 13) %>%
    filter(intervention == "Wildfire mode") %>%
    summarise(mean_red = mean(exp_red),
              sd_red = sd(exp_red)
    )

# quantify uncertainties of house size
plot_data %>%
    filter(house_size == "large") %>%
    summarise(mean_red = mean(exp_red),
              sd_red = sd(exp_red)
    )

# quantify uncertainties of house leak
prep_plot_io %>%
    select(starts_with("exp_red"), house_size, MERV, house_leak) %>%
    pivot_longer(cols = starts_with("exp_red"), names_to = "intervention", values_to = "exp_red") %>%
    filter(intervention %in% c("exp_red_shelter", "exp_red_baseline_hvac", "exp_red_baseline_out", "exp_red_baseline_pac")) %>%
    mutate(
        intervention = factor(intervention, 
                              levels = c("exp_red_shelter", "exp_red_baseline_hvac",
                                         "exp_red_baseline_out", 
                                         "exp_red_baseline_pac"),
                              labels = c("Sheltering indoors", "Air conditioning", "Wildfire mode",
                                         "Portable air cleaner")),
        house_size = factor(house_size, levels = c("small", "medium", "large")),
        MERV = factor(MERV, levels = c("8", "10", "13", "16"))
    ) %>%
    filter(house_leak == "airtight") %>%
    summarise(mean_red = mean(exp_red),
              sd_red = sd(exp_red)
    )



###### Figure A.9 New ######
# Select relevant columns and reshape data for plotting
plot_data <- prep_plot_io %>%
    select(starts_with("exp_red"), house_size, MERV) %>%
    pivot_longer(cols = starts_with("exp_red"), names_to = "intervention", values_to = "exp_red") %>%
    filter(intervention %in% c("exp_red_hvac_shelter", "exp_red_wildfire_shelter", "exp_red_pac_shelter")) %>%
    mutate(
        intervention = factor(intervention, 
                              levels = c("exp_red_hvac_shelter",
                                         "exp_red_wildfire_shelter", 
                                         "exp_red_pac_shelter"),
                              labels = c("Air conditioning", 
                                         "Wildfire mode",
                                         "Portable air cleaner")),
        house_size = factor(house_size, levels = c("small", "medium", "large")),
        MERV = factor(MERV, levels = c("8", "10", "13", "16"))
    )

# Create the boxplot
plot_data %>%
    filter(as.numeric(as.character(MERV)) %in% c(10, 13, 16)) %>%
    ggplot(., aes(x = intervention, y = exp_red, fill = interaction(MERV, house_size) )) +
    #geom_jitter(width = 0.2, size = 0.1, alpha = 0.05) +  # Add jitter for data distribution
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplots without outliers for clarity
    coord_cartesian(ylim = c(0,70)) +
    labs(
        title = "Indoor PM2.5 exposure reduction (%) relative to sheltering indoors",
        subtitle = "Scenarious grouped per house floor area and by type of filter",
        x = "Intervention Type",
        y = "Exposure Reduction (%)",
        fill = "Central air system filter type and house size"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

ggsave(filename = here("products/gfx/figureA9.png"), width = 3600, height = 1800, dpi = 300, units = "px")



###### Figure nr. air purifiers ######
prep_plot_io %>%
    filter(as.numeric(as.character(MERV)) %in% c(10, 13)) %>%
    mutate(
        nr_pacs = (pm_25_infiltration - pm_25_baseline_hvac) / (pm_25_infiltration - pm_25_pac)
    ) %>% 
    mutate(nr_pacs_ceil = ceiling(nr_pacs)) %>%
    select(MERV, house_size, house_leak, nr_pacs_ceil) %>% glimpse() %>%
    ggplot(., aes(y = nr_pacs_ceil)) +
    geom_bar()

prep_plot_io %>%
    filter(as.numeric(as.character(MERV)) %in% c(10, 13)) %>%
    mutate(
        nr_pacs_hvac = (pm_25_infiltration - pm_25_baseline_hvac) / (pm_25_infiltration - pm_25_pac),
        nr_pacs_wildfire_mode = (pm_25_infiltration - pm_25_baseline_out) / (pm_25_infiltration - pm_25_pac),
        nr_pacs_hvac_ceil = ceiling(nr_pacs_hvac),
        nr_pacs_wild_ceil = ceiling(nr_pacs_wildfire_mode)
    ) %>%
    select(MERV, house_size, house_leak, nr_pacs_hvac_ceil, nr_pacs_wild_ceil) %>%
    group_by(MERV, house_size, house_leak) %>%
    summarise(mean_nr_pacs = mean(nr_pacs_wild_ceil, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = factor(MERV), y = mean_nr_pacs, fill = house_leak)) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_wrap(~house_size) +
    labs(
        x = "MERV Rating",
        y = "Average Number of PACs (Ceiled)",
        fill = "Leakage",
        title = "Grouped Bar Plot of PACs by MERV, House Leakiness, and Size"
    ) +
    theme_minimal()   

prep_plot_io %>%
    filter(as.numeric(as.character(MERV)) %in% c(10, 13)) %>%
    mutate(
        nr_pacs_hvac = (pm_25_infiltration - pm_25_baseline_hvac) / (pm_25_infiltration - pm_25_pac),
        nr_pacs_wildfire_mode = (pm_25_infiltration - pm_25_baseline_out) / (pm_25_infiltration - pm_25_pac),
        nr_pacs_hvac_ceil = ceiling(nr_pacs_hvac),
        nr_pacs_wild_ceil = ceiling(nr_pacs_wildfire_mode),
        house_id = paste(house_size, house_leak, sep = "_")
    ) %>%
    select(MERV, house_id, nr_pacs_hvac_ceil, nr_pacs_wild_ceil) %>%
    pivot_longer(cols = c(nr_pacs_hvac_ceil, nr_pacs_wild_ceil),
                 names_to = "mode", values_to = "nr_pacs") %>%
    group_by(MERV, house_id, mode) %>%
    summarise(mean_nr_pacs = mean(nr_pacs, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = interaction(factor(MERV), house_id), y = mean_nr_pacs, fill = mode)) +
    geom_col(position = position_dodge(width = 0.8)) +
    labs(
        x = "MERV + House (Size_Leak)",
        y = "Average Number of PACs (Ceiled)",
        fill = "Operation Mode",
        title = "Comparison of PAC Needs by MERV, House, and Operation Mode"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### others #####








df_sim_plt %>%
    #filter(MERV %in% c(8,13)) %>%
    group_by(test) %>%
    summarise(pm_25_infiltration = sum(pm_25_infiltration),
              fan_10min_hvac = sum(fan_10min_hvac),
              pm_25_10min_hvac = sum(pm_25_10min_hvac),
              fan_baseline_hvac = sum(fan_baseline_hvac),
              pm_25_baseline_hvac = sum(pm_25_baseline_hvac),
              fan_baseline_out = sum(fan_baseline_out),
              pm_25_baseline_out = sum(pm_25_baseline_out),
              fan_baseline_in = sum(fan_baseline_in),
              pm_25_baseline_in = sum(pm_25_baseline_in),
              fan_AQ_logic = sum(fan_AQ_logic),
              pm_25_AQ_logic = sum(pm_25_AQ_logic),
              fan_pac = sum(fan_pac),
              pm_25_pac = sum(pm_25_pac),
              pm_25_out = sum(pm_25_out),
              geoid = first(geoid),
              start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage)
    ) %>% 
    mutate(exp_red_shelter = 100 - pm_25_infiltration / pm_25_out * 100,
           exp_red_10_min = 100 - pm_25_10min_hvac / pm_25_out * 100,
           exp_red_baseline_hvac = 100 - pm_25_baseline_hvac / pm_25_out * 100,
           exp_red_baseline_out = 100 - pm_25_baseline_out / pm_25_out * 100,
           exp_red_baseline_in = 100 - pm_25_baseline_in / pm_25_out * 100,
           exp_red_baseline_AQ = 100 - pm_25_AQ_logic / pm_25_out * 100,
           exp_red_baseline_pac = 100 - pm_25_pac / pm_25_out * 100
    ) %>% glimpse() %>%
    ungroup() %>%
    mutate(house_size = case_when(house_area <= small_house ~ "small",
                                  house_area > small_house & house_area <= large_house ~ "medium",
                                  house_area > large_house ~ "large"),
           house_leak = case_when(effect_leakage <= airtight_house ~ "airtight",
                                  effect_leakage > airtight_house & effect_leakage <= leaky_house ~ "med. leaky",
                                  effect_leakage > leaky_house ~ "very leaky"
           ),
           # Combine factors and ensure proper ordering
           house_combined = interaction(house_size, house_leak),
           house_size = factor(house_size, levels = c("small", "medium", "large")), # Ensures correct order
           house_leak = factor(house_leak, levels = c("airtight", "med. leaky", "very leaky")), # Ensures correct order
           house_combined = factor(house_combined, levels = c(
               "small.airtight", "small.med. leaky", "small.very leaky",
               "medium.airtight", "medium.med. leaky", "medium.very leaky",
               "large.airtight", "large.med. leaky", "large.very leaky"
           ))
    ) %>% glimpse() %>%
    select(test, geoid, start_IO, MERV, house_area, house_size, effect_leakage, house_leak, 
           house_combined,
           exp_red_shelter, 
           #exp_red_10_min, 
           exp_red_baseline_hvac,
           exp_red_baseline_out,
           #exp_red_baseline_in,
           #exp_red_baseline_AQ, 
           exp_red_baseline_pac
    ) %>% glimpse() %>%
    pivot_longer(
        cols = starts_with("exp_red"),
        names_to = "exp_red_type",
        values_to = "exp_red_value"
    ) %>% glimpse() %>%
    #ggplot(., aes(x = house_size, y = exp_red_value, fill = factor(exp_red_type), color = factor(MERV))) +
    # Boxplot with Ordered House Size and House Leak
    #. ggplot(., aes(x = house_combined, y = exp_red_value, fill = factor(MERV))) +
    #. geom_boxplot() +
    #. facet_wrap(~ factor(exp_red_type)) +
    #. labs(title = "Boxplot with Ordered House Size and House Leak",
    #.      x = "House Size & Leakiness Combination",
    #.      y = "Exposure Reduction Value",
    #.      fill = "MERV Rating") +
    #. theme_minimal() +
    #. theme(axis.text.x = element_text(angle = 45, hjust = 1),
#.       panel.background = element_rect(fill = "white", color = "black"),  # White panel with black border
#.       plot.background = element_rect(fill = "white", color = NA)  # White background, no border
#.       )
# Boxplot with House Size and House Leak Facets
ggplot(., aes(x = factor(house_size), y = exp_red_value, fill = factor(exp_red_type))) +
    geom_boxplot() +
    facet_wrap(~ factor(house_leak)) +
    labs(title = "Boxplot with House Size and House Leak Facets",
         x = "House Size",
         y = "Exposure Reduction Value",
         fill = "Intervention Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", color = "black"),  # White panel with black border
          plot.background = element_rect(fill = "white", color = NA)  # White background, no border
    )

#ggsave(filename = here("products/gfx/sens_analysis_comb_size_leak.png"), width = 3600, height = 1800, dpi = 300, units = "px")

ggsave(filename = here("products/gfx/sens_analysis_size_leak.png"), width = 3600, height = 1800, dpi = 300, units = "px")



# Boxplot with House Size and House Leak Facets

# Define custom labels for facet titles
exp_red_labels <- c(
    airtight = "Low Level of Air Infiltration",
    `med. leaky` = "Moderate Level of Air Infiltration",
    `very leaky` = "High Level of Air Infiltration"
)

df_sim_plt %>%
    #filter(MERV %in% c(8,13)) %>%
    group_by(test) %>%
    summarise(pm_25_infiltration = sum(pm_25_infiltration),
              fan_10min_hvac = sum(fan_10min_hvac),
              pm_25_10min_hvac = sum(pm_25_10min_hvac),
              fan_baseline_hvac = sum(fan_baseline_hvac),
              pm_25_baseline_hvac = sum(pm_25_baseline_hvac),
              fan_baseline_out = sum(fan_baseline_out),
              pm_25_baseline_out = sum(pm_25_baseline_out),
              fan_baseline_in = sum(fan_baseline_in),
              pm_25_baseline_in = sum(pm_25_baseline_in),
              fan_AQ_logic = sum(fan_AQ_logic),
              pm_25_AQ_logic = sum(pm_25_AQ_logic),
              fan_pac = sum(fan_pac),
              pm_25_pac = sum(pm_25_pac),
              pm_25_out = sum(pm_25_out),
              geoid = first(geoid),
              start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage)
    ) %>% 
    mutate(exp_red_shelter = 100 - pm_25_infiltration / pm_25_out * 100,
           exp_red_10_min = 100 - pm_25_10min_hvac / pm_25_out * 100,
           exp_red_baseline_hvac = 100 - pm_25_baseline_hvac / pm_25_out * 100,
           exp_red_baseline_out = 100 - pm_25_baseline_out / pm_25_out * 100,
           exp_red_baseline_in = 100 - pm_25_baseline_in / pm_25_out * 100,
           exp_red_baseline_AQ = 100 - pm_25_AQ_logic / pm_25_out * 100,
           exp_red_baseline_pac = 100 - pm_25_pac / pm_25_out * 100
    ) %>% glimpse() %>%
    ungroup() %>%
    mutate(house_size = case_when(house_area <= small_house ~ "small",
                                  house_area > small_house & house_area <= large_house ~ "medium",
                                  house_area > large_house ~ "large"),
           house_leak = case_when(effect_leakage <= airtight_house ~ "airtight",
                                  effect_leakage > airtight_house & effect_leakage <= leaky_house ~ "med. leaky",
                                  effect_leakage > leaky_house ~ "very leaky"
           ),
           # Combine factors and ensure proper ordering
           house_combined = interaction(house_size, house_leak),
           house_size = factor(house_size, levels = c("small", "medium", "large")), # Ensures correct order
           house_leak = factor(house_leak, levels = c("airtight", "med. leaky", "very leaky")), # Ensures correct order
           house_combined = factor(house_combined, levels = c(
               "small.airtight", "small.med. leaky", "small.very leaky",
               "medium.airtight", "medium.med. leaky", "medium.very leaky",
               "large.airtight", "large.med. leaky", "large.very leaky"
           ))
    ) %>% glimpse() %>%
    select(test, geoid, start_IO, MERV, house_area, house_size, effect_leakage, house_leak, 
           house_combined,
           exp_red_shelter, 
           #exp_red_10_min, 
           exp_red_baseline_hvac,
           exp_red_baseline_out,
           #exp_red_baseline_in,
           #exp_red_baseline_AQ, 
           exp_red_baseline_pac
    ) %>% glimpse() %>%
    pivot_longer(
        cols = starts_with("exp_red"),
        names_to = "exp_red_type",
        values_to = "exp_red_value"
    ) %>% glimpse() %>%
    ggplot(., aes(x = factor(house_size), y = exp_red_value, fill = factor(exp_red_type))) +
    geom_boxplot() +
    facet_wrap(~ factor(house_leak), labeller = as_labeller(exp_red_labels)) +
    labs(title = "Boxplot with House Size and House Leak Facets",
         x = "House Size",
         y = "Exposure Reduction Value",
         fill = "Intervention Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", color = "black"),  # White panel with black border
          plot.background = element_rect(fill = "white", color = NA)  # White background, no border
    )

ggsave(filename = here("products/gfx/sens_analysis_size_leak.png"), width = 3600, height = 1800, dpi = 300, units = "px")


##### Boxplot with Ordered House Size and House Leak #####

# Define custom labels for facet titles
exp_red_labels <- c(
    exp_red_shelter = "Sheltering Indoors Exposure Reduction (%)",
    exp_red_baseline_hvac = "Air Conditioning Exposure Reduction (%)",
    exp_red_baseline_out = "Wildfire Mode Exposure Reduction (%)",
    exp_red_baseline_pac = "Portable Air Cleaner Exposure Reduction (%)"
)

df_sim_plt %>%
    #filter(MERV %in% c(8,13)) %>%
    group_by(test) %>%
    summarise(pm_25_infiltration = sum(pm_25_infiltration),
              fan_10min_hvac = sum(fan_10min_hvac),
              pm_25_10min_hvac = sum(pm_25_10min_hvac),
              fan_baseline_hvac = sum(fan_baseline_hvac),
              pm_25_baseline_hvac = sum(pm_25_baseline_hvac),
              fan_baseline_out = sum(fan_baseline_out),
              pm_25_baseline_out = sum(pm_25_baseline_out),
              fan_baseline_in = sum(fan_baseline_in),
              pm_25_baseline_in = sum(pm_25_baseline_in),
              fan_AQ_logic = sum(fan_AQ_logic),
              pm_25_AQ_logic = sum(pm_25_AQ_logic),
              fan_pac = sum(fan_pac),
              pm_25_pac = sum(pm_25_pac),
              pm_25_out = sum(pm_25_out),
              geoid = first(geoid),
              start_IO = first(start_IO),
              MERV = first(MERV),
              house_area = first(house_area),
              effect_leakage = first(effect_leakage)
    ) %>% 
    mutate(exp_red_shelter = 100 - pm_25_infiltration / pm_25_out * 100,
           exp_red_10_min = 100 - pm_25_10min_hvac / pm_25_out * 100,
           exp_red_baseline_hvac = 100 - pm_25_baseline_hvac / pm_25_out * 100,
           exp_red_baseline_out = 100 - pm_25_baseline_out / pm_25_out * 100,
           exp_red_baseline_in = 100 - pm_25_baseline_in / pm_25_out * 100,
           exp_red_baseline_AQ = 100 - pm_25_AQ_logic / pm_25_out * 100,
           exp_red_baseline_pac = 100 - pm_25_pac / pm_25_out * 100
    ) %>% glimpse() %>%
    ungroup() %>%
    mutate(house_size = case_when(house_area <= small_house ~ "small",
                                  house_area > small_house & house_area <= large_house ~ "medium",
                                  house_area > large_house ~ "large"),
           house_leak = case_when(effect_leakage <= airtight_house ~ "airtight",
                                  effect_leakage > airtight_house & effect_leakage <= leaky_house ~ "med. leaky",
                                  effect_leakage > leaky_house ~ "very leaky"
           ),
           # Combine factors and ensure proper ordering
           house_combined = interaction(house_size, house_leak),
           house_size = factor(house_size, levels = c("small", "medium", "large")), # Ensures correct order
           house_leak = factor(house_leak, levels = c("airtight", "med. leaky", "very leaky")), # Ensures correct order
           house_combined = factor(house_combined, levels = c(
               "small.airtight", "small.med. leaky", "small.very leaky",
               "medium.airtight", "medium.med. leaky", "medium.very leaky",
               "large.airtight", "large.med. leaky", "large.very leaky"
           ))
    ) %>% glimpse() %>%
    select(test, geoid, start_IO, MERV, house_area, house_size, effect_leakage, house_leak, 
           house_combined,
           exp_red_shelter, 
           #exp_red_10_min, 
           exp_red_baseline_hvac,
           exp_red_baseline_out,
           #exp_red_baseline_in,
           #exp_red_baseline_AQ, 
           exp_red_baseline_pac
    ) %>% glimpse() %>%
    pivot_longer(
        cols = starts_with("exp_red"),
        names_to = "exp_red_type",
        values_to = "exp_red_value"
    ) %>% glimpse() %>%
    ggplot(., aes(x = house_combined, y = exp_red_value, fill = factor(MERV))) +
    geom_boxplot() +
    facet_wrap(~ factor(exp_red_type), labeller = as_labeller(exp_red_labels)) +
    labs(title = "Boxplot with Ordered House Size and House Leak",
         x = "House Size & Leakiness Combination",
         y = "Exposure Reduction Value",
         fill = "MERV Rating") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", color = "black"),  # White panel with black border
          plot.background = element_rect(fill = "white", color = NA)  # White background, no border
    )

ggsave(filename = here("products/gfx/sens_analysis_comb_size_leak.png"), width = 3600, height = 1800, dpi = 300, units = "px")