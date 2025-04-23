# Using smart thermostats to reduce indoor exposure to wildfire fine particulate matter (PM2.5)
# Cite: https://doi.org/10.1016/j.indenv.2025.100088
# ----
# TASK: Dynamic PM2.5 by Dynamic Infiltration and County
# Code Authors: Federico Dallo, Thomas Parkinson, Carlos Duarte, Chai Yoon Um, Paul Raftery
# ----

#### LIBRARY SETUP ####

library(pacman)
pacman::p_load(
    tidyverse, here, lubridate, zoo,
    pbapply, furrr, future.apply
)

here::i_am("tstat-wildfire.Rproj")


#### DATA LOADING ####

# MERV filter efficiency and pressure loss
df_filter <- tibble(
    filter_merv = c(0, 8, 9, 10, 11, 12, 13, 14, 15, 16),
    filter_eff  = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 0.85, 0.9, 0.9, 0.95),
    filter_pres = c(0, 0.056, 0.056, 0.056, 0.088, 0.088, 0.088, 0.120, 0.120, 0.120)
)

# Load tract, ISD, EPA, Ecobee data
# Note by F. Dallo: data will be soon provided in a separate repo. Request data via email if needed: federico.dallo@cnr.it
df_tract <- read_rds(here("data", "US", "CA", "df_tract.rds"))
df_isd   <- read_rds(here("data", "US", "CA", "df_isd.rds"))
df_epa   <- read_rds(here("data", "US", "CA", "df_epa.rds"))
df_ecobee <- read_rds(here("data", "US", "CA", "df_ecobee.rds"))
df_counties_baseline_run <- read_rds(here("data", "US", "CA", "df_counties_baseline_run.rds"))

# Expand time and interpolate for ISD and EPA
df_isd <- df_isd %>%
    group_by(isd_station) %>%
    arrange(datetime) %>%
    complete(datetime = seq(min(datetime), max(datetime), by = "10 min")) %>%
    mutate(
        wind_speed = na.approx(wind_speed, rule = 2),
        t_out = na.approx(t_out, rule = 2)
    ) %>%
    ungroup()

df_epa <- df_epa %>%
    mutate(pm_25_out = ifelse(is.nan(pm_25_out), NA, pm_25_out)) %>%
    group_by(epa_station) %>%
    arrange(datetime) %>%
    complete(datetime = seq(min(datetime), max(datetime), by = "10 min")) %>%
    mutate(pm_25_out = na.approx(pm_25_out, rule = 2)) %>%
    ungroup()

df_ecobee <- df_ecobee %>%
    drop_na(county) %>%
    group_by(date_time, county) %>%
    summarise(t_in = mean(t_ctrl, na.rm = TRUE), .groups = "drop") %>%
    group_by(county) %>%
    arrange(date_time) %>%
    complete(date_time = seq(min(date_time), max(date_time), by = "10 min")) %>%
    mutate(t_in = na.approx(t_in, rule = 2)) %>%
    ungroup()


#### TRACT PREPARATION ####

# Assign heat load based on ASHRAE climate zone
df_tract <- df_tract %>%
    mutate(heat_load = case_when(
        climate_ashrae == "2B" ~ 40,
        climate_ashrae == "3B" ~ 45,
        climate_ashrae == "3C" ~ 45,
        climate_ashrae == "4B" ~ 50,
        climate_ashrae == "5B" ~ 60,
        climate_ashrae == "6B" ~ 65
    ))

# Map counties without indoor temp data to nearby ones
df_tract <- df_tract %>%
    mutate(county = case_when(
        county %in% c("colusa", "glenn", "tehama") ~ "sutter",
        county %in% c("del norte", "mendocino", "trinity") ~ "humboldt",
        county == "mariposa" ~ "tuolumne",
        county == "modoc" ~ "lassen",
        TRUE ~ county
    ))

# Compute building parameters
df_tract <- df_tract %>%
    mutate(
        house_volume = house_area * house_height,
        floor_area_ft2 = house_area * 10.764,
        pred_heat_load_btu_hr = round((floor_area_ft2 * heat_load) / 3000) * 3000,
        pred_air_flow_rate_cfm_hi = round((pred_heat_load_btu_hr / (60 * 0.0745 * 0.24025 * 70)) / 50) * 50,
        design_furnace = pred_air_flow_rate_cfm_hi * (0.3048^3)
    ) %>%
    select(geoid, county, epa_station, isd_station, house_height, house_volume, house_area, effect_leakage, design_furnace)


#### DYNAMIC MODEL CONSTANTS ####

model_const <- lst(
    c_g = 9.8,
    c_r = 0.5,
    c_t0 = 298,
    c_x = 0.25,
    c_c = 0.19,
    c_factor = 10,
    c_a = 0.67,
    c_b = 0.25
)

if (!dir.exists(here("data", "simulations"))) {
    dir.create(here("data", "simulations"))
}

list_sim_file <- list.files(here("data", "simulations"), pattern = "[0-9].rds") %>%
    str_replace(".rds", "")

df_tract_reduced <- df_tract %>%
    filter(!geoid %in% list_sim_file)


#### PARALLEL SIMULATION ####

plan(multisession, workers = parallel::detectCores() - 1)

results <- future_map(seq_len(nrow(df_tract_reduced)), process_tract, .progress = TRUE)

stop()


#### FUNCTION: process_tract() ####

process_tract <- function(i) { # for parallelization
    
    print(nrow(df_tract_reduced) + 1 - i)
    
    # get tract info
    sim_house <- df_tract_reduced %>%
        slice(i)
    
    # print 
    print(sim_house$geoid)
    
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
        select(geoid, county, datetime, d_t, house_volume, house_area, v_infil, design_furnace, pm_25_out)
    
    # add the ecobee county fan runtime
    sim_house <- sim_house %>%
        left_join(., df_counties_baseline_run, by = c("county", "datetime" = "date_time"))
    
    # define model variables 
    sim_house <- sim_house %>%
        mutate(pm25_pene = 0.75, # penetration factor of PM 2.5 particles
               pm25_dep = 0.5, # [- per hour] deposition per hour
               filter_merv = 13, # merv rating
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
            c_in_baseline = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # estimate initial indoor pm2.5 concentration (#FEDE 0.5 instead of 0.35 - be conservative)
            v_fan_baseline = ifelse(minute(datetime) < 0, yes = 1, no = 0) * design_furnace * d_t, # m3
            v_da_baseline = v_fan_baseline * ((1 - filter_pres) - n_ls),
            v_ra_baseline = (v_infil + v_da_baseline) * as.integer(v_fan_baseline > 0), # m3
            c_bf_baseline = 0, # (v_ra_baseline*c_in_baseline + v_fan_baseline*n_lr*pm_25_out) / (v_ra_baseline + v_fan_baseline*n_lr), 
            c_af_baseline = c_bf_baseline*(1 - filter_eff),
            
            # 10-min hvac 
            c_in_hvac = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration
            v_fan_hvac = ifelse(minute(datetime) < 10, yes = 1, no = 0) * design_furnace * d_t, # m3 ##FEDE remove hard coded 10.. 
            v_da_hvac = v_fan_hvac * ((1 - filter_pres) - n_ls), # m3
            v_ra_hvac = (v_infil + v_da_hvac) * as.integer(v_fan_hvac > 0), # m3, second part evaluating that is non-zero only when HVAC is ON
            c_bf_hvac = case_when( (v_ra_hvac + v_fan_hvac*n_lr) > 0 ~ (v_ra_hvac*c_in_hvac + v_fan_hvac*n_lr*pm_25_out)/(v_ra_hvac + v_fan_hvac*n_lr), # ug/m3
                                   (v_ra_hvac + v_fan_hvac*n_lr) == 0 ~ 0,
                                   TRUE ~ 0
            ),
            c_af_hvac = c_bf_hvac*(1 - filter_eff), # ug/m3
            
            # baseline county hvac NEW
            c_in_base_hvac = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
            v_fan_base_hvac = fan_baseline_logic * design_furnace * d_t, # m3 
            v_da_base_hvac = v_fan_base_hvac * ((1 - filter_pres) - n_ls), # m3
            v_ra_base_hvac = (v_infil + v_da_base_hvac) * as.integer(v_fan_base_hvac > 0), # m3, second part evaluating that is non-zero only when HVAC is ON
            c_bf_base_hvac = case_when( (v_ra_base_hvac + v_fan_base_hvac*n_lr) > 0 ~ (v_ra_base_hvac*c_in_base_hvac + v_fan_base_hvac*n_lr*pm_25_out)/(v_ra_base_hvac + v_fan_base_hvac*n_lr), # ug/m3
                                        (v_ra_base_hvac + v_fan_base_hvac*n_lr) == 0 ~ 0,
                                        TRUE ~ 0
            ),
            c_af_base_hvac = c_bf_base_hvac*(1 - filter_eff), # ug/m3
            
            # control logic hvac with BASELINE and outdoor PM2.5 >= 35
            c_in_base_out_logic = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
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
            c_in_base_in_logic = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
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
            c_in_logic = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # ug/m3 estimate initial indoor pm2.5 concentration 
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
            c_in_pac = if_else(row_number() == 1, pm_25_out*0.5, NA_real_), # estimate initial indoor pm2.5 concentration 
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
            
            
            # HVAC control logic over baseline and when outdoor PM2.5 >= 35 ug/m3 TODO
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
            
            # HVAC control logic over baseline and when indoor PM2.5 >= 5 ug/m3 TODO
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
        select(geoid, datetime, pm_25_out,
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
    write_rds(sim_house, file = here("data", "simulations", paste0(sim_house[1,1],'.rds')))
    
    # summarise exposure and write to file # FEDE check runtime
    sim_house %>%
        summarise(geoid = first(geoid),
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
        write_delim(., here("data", "simulations", "_total_exposure.csv"), delim = ",", append = TRUE)
    
    return(sim_house)
}


#### IMPORT AND SUMMARIZE OUTPUT ####

list_file <- list.files(path = here("data", "simulations"), pattern = "[0-9].rds", full.names = TRUE)

import_sim <- function(sim_file) {
    read_rds(file = sim_file) %>%
        group_by(datetime = floor_date(datetime, unit = "hour")) %>%
        summarise(
            geoid = first(geoid),
            pm_25_infiltration = mean(pm_25_infiltration, na.rm = TRUE),
            fan_10min_hvac = mean(fan_10min_hvac, na.rm = TRUE) * 10,
            pm_25_10min_hvac = mean(pm_25_10min_hvac, na.rm = TRUE),
            fan_baseline_hvac = mean(fan_baseline_hvac, na.rm = TRUE) * 10,
            pm_25_baseline_hvac = mean(pm_25_baseline_hvac, na.rm = TRUE),
            fan_baseline_out = mean(fan_baseline_out, na.rm = TRUE) * 10,
            pm_25_baseline_out = mean(pm_25_baseline_out, na.rm = TRUE),
            fan_baseline_in = mean(fan_baseline_in, na.rm = TRUE) * 10,
            pm_25_baseline_in = mean(pm_25_baseline_in, na.rm = TRUE),
            fan_AQ_logic = mean(fan_AQ_logic, na.rm = TRUE) * 10,
            pm_25_AQ_logic = mean(pm_25_AQ_logic, na.rm = TRUE),
            fan_pac = mean(fan_pac, na.rm = TRUE) * 10,
            pm_25_pac = mean(pm_25_pac, na.rm = TRUE),
            pm_25_out = mean(pm_25_out, na.rm = TRUE)
        )
}

list_sim <- pblapply(list_file, import_sim)
df_sim <- bind_rows(list_sim)

write_rds(df_sim, here("data", "simulations", paste0("df_hourly_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")), compress = "gz")


#### CLEANUP ####

rm(list_sim, list_file, import_sim)