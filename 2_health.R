# Using smart thermostats to reduce indoor exposure to wildfire fine particulate matter (PM2.5)
# Cite: https://doi.org/10.1016/j.indenv.2025.100088
# ----
# TASK: Health Impact Assessment
# Code Authors: Federico Dallo, Thomas Parkinson, Carlos Duarte, Chai Yoon Um, Paul Raftery
# ----


#### SETUP ####

# use pacman for package management
library(pacman)

# load packages
pacman::p_load(tidyverse, here, lubridate, # essential
               scales, patchwork, # plots
               usmap, sf, tigris, # mapping
               #PWFSLSmoke,
               tidygeocoder, # EPA AirNow pm2.5 data
               gganimate) # animate map

# Use "here" package to set working directory
here::i_am("healthRiskADAPT.Rproj")

# turn off scientific notation
options(scipen = 999)



#### LOAD SIMULATION DATA ####

# load tract dataset
df_tract <- read_rds(here("data", "US", "CA", "df_tract.rds"))

# load simulated exposure summary file and add to dataset
df_tract <- read_csv(here("data", "simulations", "_total_exposure.csv"), 
                     col_names = c("geoid",
                                   "samples",
                                   "pm_25_infiltration",
                                   "fan_10min_hvac",
                                   "pm_25_10min_hvac",
                                   "fan_baseline_hvac",
                                   "pm_25_baseline_hvac",
                                   "fan_baseline_out",
                                   "pm_25_baseline_out",
                                   "fan_baseline_in",
                                   "pm_25_baseline_in",
                                   "fan_AQ_logic",
                                   "pm_25_AQ_logic",
                                   "fan_pac",
                                   "pm_25_pac",
                                   "pm_25_out"),
                     col_types = "cdddddd") %>%
    mutate(pm_25_diff_baseline = (pm_25_baseline_hvac - pm_25_infiltration) / (samples),
           pm_25_diff_logic = (pm_25_baseline_out - pm_25_infiltration) / (samples)) %>%
    distinct() %>%
    left_join(df_tract, ., by = "geoid")




#### MORTALITY ####

# Health impact based on benmap and Zhao et al. 
## https://www.epa.gov/sites/default/files/2015-04/documents/benmap-ce_user_manual_march_2015.pdf
## https://www.mdpi.com/1660-4601/12/7/8448/htm#B7-ijerph-12-08448

# health impact assessment
# All cause mortality, short-term exposure, daily (24h), log-linear page 291 in BenMAP manual
## average_pm = pm_25_in * 0.9 + pm_25_out * 0.1 
## incidence_daily = (1-(1/exp(-1 * beta * deltaPM))) * incidence * population * convert_annual_to_daily #log-linear function from BenMAP, EPA
## beta = 0.000145 #Ito, 2013. page 371 in BenMAP manual
## standard_error = 0.000075 #Ito, 2013. page 371 in BenMAP manual
## beta_lower = beta - 1.96 * standard_error
## beta_upper = beta + 1.96 * standard_error
## deltaPM = average_pm - 0
## incidence = 7.4 / 1000 # https://ehp.niehs.nih.gov/doi/full/10.1289/ehp.1104035
## convert_annual_to_daily = 1/365

# load hourly dataset and summarise daily values
list_file_df_hourly <- list.files(path = here("data", "simulations"), pattern = "df_hourly", full.names = FALSE)

recent_df_hourly <- list_file_df_hourly %>% 
    tail(1)

df_daily <- read_rds(here("data", "simulations", recent_df_hourly)) %>%
    group_by(geoid, 
             date = floor_date(datetime, "1 day")) %>%
    summarise(samples = n(),
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
              pm_25_out = sum(pm_25_out, na.rm = TRUE)) %>%
    ungroup()

# add in population
df_daily <- df_tract %>%
    select(geoid, population) %>%
    left_join(df_daily, ., by = "geoid")

# calculate average exposed concentration (90% inside, 10% outside)
df_health <- df_daily %>%
    mutate(pm_25_infiltration = pm_25_infiltration / samples,
           pm_25_10min_hvac = pm_25_10min_hvac / samples,
           pm_25_baseline_hvac = pm_25_baseline_hvac / samples,
           pm_25_baseline_out = pm_25_baseline_out / samples,
           pm_25_baseline_in = pm_25_baseline_in / samples,
           pm_25_AQ_logic = pm_25_AQ_logic / samples,
           pm_25_pac = pm_25_pac / samples,
           pm_25_out = pm_25_out / samples) %>%
    mutate(pm_25_infiltration = pm_25_infiltration * 0.9 + pm_25_out * 0.1,
           pm_25_10min_hvac = pm_25_10min_hvac * 0.9 + pm_25_out * 0.1,
           pm_25_baseline_hvac = pm_25_baseline_hvac * 0.9 + pm_25_out * 0.1,
           pm_25_baseline_out = pm_25_baseline_out * 0.9 + pm_25_out * 0.1,
           pm_25_baseline_in = pm_25_baseline_in * 0.9 + pm_25_out * 0.1,
           pm_25_AQ_logic = pm_25_AQ_logic * 0.9 + pm_25_out * 0.1,
           pm_25_pac = pm_25_pac * 0.9 + pm_25_out * 0.1
    )

# calculate mortality risk
df_health <- df_health %>%
    mutate(cut_off = 0,
           beta = 0.000145,
           standard_error = 0.000075,
           beta_lower = beta - 1.96 * standard_error,
           beta_upper = beta + 1.96 * standard_error,
           incidence = 7.38 * 10^(-3),
           convert_annual_to_daily = 1/365,
           days = samples / 24,
           pm0 = 10, # WHO value
           # infiltration scenario
           mortality_infiltration = (1 - exp(beta * (pm0 - pm_25_infiltration))) * incidence * convert_annual_to_daily,
           incidence_infiltration = ifelse(mortality_infiltration * population <= 0, 0, mortality_infiltration * population),
           # 10-min logic scenario
           mortality_10min_hvac = (1 - exp(beta * (pm0 - pm_25_10min_hvac))) * incidence * convert_annual_to_daily,
           incidence_10min_hvac = ifelse(mortality_10min_hvac * population <= 0, 0, mortality_10min_hvac * population),
           # baseline hvac logic scenario
           mortality_baseline_hvac = (1 - exp(beta * (pm0 - pm_25_baseline_hvac))) * incidence * convert_annual_to_daily,
           incidence_baseline_hvac = ifelse(mortality_baseline_hvac * population <= 0, 0, mortality_baseline_hvac * population),
           # baseline outdoor PM2.5 logic scenario
           mortality_baseline_out = (1 - exp(beta * (pm0 - pm_25_baseline_out))) * incidence * convert_annual_to_daily,
           incidence_baseline_out = ifelse(mortality_baseline_out * population <= 0, 0, mortality_baseline_out * population),
           # baseline indoor PM2.5 logic scenario
           mortality_baseline_in = (1 - exp(beta * (pm0 - pm_25_baseline_in))) * incidence * convert_annual_to_daily,
           incidence_baseline_in = ifelse(mortality_baseline_in * population <= 0, 0, mortality_baseline_in * population),
           # AQ logic scenario
           mortality_AQ = (1 - exp(beta * (pm0 - pm_25_AQ_logic))) * incidence * convert_annual_to_daily,
           incidence_AQ = ifelse(mortality_AQ * population <= 0, 0, mortality_AQ * population),
           # PAC logic scenario
           mortality_pac = (1 - exp(beta * (pm0 - pm_25_pac))) * incidence * convert_annual_to_daily,
           incidence_pac = ifelse(mortality_pac * population <= 0, 0, mortality_pac * population)
    )

# calculate % reduction in excess mortality from HVAC and value of mortality risk reduction
df_perc_red <- df_health %>%
    group_by(geoid) %>%
    summarise(# infiltration scenario
        infiltration_scenario = sum(incidence_infiltration, na.rm = TRUE),
        # 10-min logic scenario
        fix_min_hvac_scenario = sum(incidence_10min_hvac, na.rm = TRUE),
        # baseline hvac logic scenario
        baseline_hvac_scenario = sum(incidence_baseline_hvac, na.rm = TRUE),
        # baseline outdoor PM2.5 logic scenario
        baseline_out_scenario = sum(incidence_baseline_out, na.rm = TRUE),
        # baseline indoor PM2.5 logic scenario
        baseline_in_scenario = sum(incidence_baseline_in, na.rm = TRUE),
        # AQ logic scenario
        AQ_scenario = sum(incidence_AQ, na.rm = TRUE),
        # PAC logic scenario
        pac_scenario = sum(incidence_pac, na.rm = TRUE),
        population = median(population, na.rm = TRUE)) %>%
    #### old
    #baseline_scenario = sum(incidence_baseline, na.rm = TRUE),
    #hvac_scenario = sum(incidence_hvac, na.rm = TRUE),
    #logic_scenario = sum(incidence_logic, na.rm = TRUE),
    #population = median(population, na.rm = TRUE)) %>%
    ungroup() %>%
    summarise(infiltration_scenario = sum(infiltration_scenario, na.rm = TRUE),
              fix_min_hvac_scenario = sum(fix_min_hvac_scenario, na.rm = TRUE),
              baseline_hvac_scenario = sum(baseline_hvac_scenario, na.rm = TRUE),
              baseline_out_scenario = sum(baseline_out_scenario, na.rm = TRUE),
              baseline_in_scenario = sum(baseline_in_scenario, na.rm = TRUE),
              AQ_scenario = sum(AQ_scenario, na.rm = TRUE),
              pac_scenario = sum(pac_scenario, na.rm = TRUE),
              population = sum(population, na.rm = TRUE)) %>%
    mutate(perc_reduction_baseline = scales::percent(1 - (baseline_hvac_scenario/infiltration_scenario)),
           perc_reduction_fixmin_hvac = scales::percent(1 - (fix_min_hvac_scenario/baseline_hvac_scenario)),
           perc_reduction_baseline_out = scales::percent(1 - (baseline_out_scenario/baseline_hvac_scenario)),
           perc_reduction_baseline_in = scales::percent(1 - (baseline_in_scenario/baseline_hvac_scenario)),
           perc_reduction_AQ = scales::percent(1 - (AQ_scenario/baseline_hvac_scenario)),
           perc_reduction_pac = scales::percent(1 - (pac_scenario/baseline_hvac_scenario)),
           vsl_reduction_logic = (baseline_hvac_scenario - baseline_out_scenario) * 9960000, # based on 7.2M in 2006 adjusted for inflation
    )
# old
#perc_reduction_hvac = scales::percent(1 - (hvac_scenario/baseline_scenario)),
#perc_reduction_logic = scales::percent(1 - (logic_scenario/baseline_scenario)),
#vsl_reduction_logic = (baseline_scenario - hvac_scenario) * 9960000, # based on 7.2M in 2006 adjusted for inflation
#vsl_per_person = vsl_reduction_logic / population)

# clean up
rm(df_health, df_daily)


#### VISUALIZATION ####

## MAP ##

# load county dataset
df_county <- tigris::counties(state = "California", cb = TRUE) %>%
    st_as_sf()

# make list of norcal tracts
list_tracts <- df_county %>%
    mutate(name = str_to_lower(NAME)) %>%
    filter(!name %in% c("san luis obispo", "kern", "san bernardino", "santa barbara", "ventura", 
                        "los angeles", "orange", "riverside", "san diego", "imperial"))

# get wildfire data from https://data-nifc.opendata.arcgis.com/datasets/nifc::wfigs-wildland-fire-locations-full-history/explore
df_fires <- read_rds(here("data", "US", "df_fires.rds")) %>%
    filter(state == "US-CA",
           start >= "2020-08-01" & start < "2020-10-01",
           size > 100000,
           str_to_lower(county) %in% list_tracts$name)


# map indoor PM2.5 concentrations by tract
p1 <- ggplot() + 
    geom_sf(data = st_as_sf(df_tract), aes(fill = pm_25_diff_logic), alpha = 0.90, colour = NA, size = 0.05) +
    geom_sf(data = df_county, fill = NA, colour = "white", size = 0.05) +
    geom_point(data = df_fires, aes(x = long, y = lat), colour = "#FBE40E", size = 6, alpha = 0.25, shape = 16) +
    geom_point(data = df_fires, aes(x = long, y = lat), colour = "#FBE40E", size = 2.5, alpha = 0.5, shape = 16) +
    geom_point(data = df_fires, aes(x = long, y = lat), colour = "#FBE40E", size = 1, alpha = 0.8, shape = 16) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_size_continuous(limits = c(0, 4)) +
    scale_fill_viridis_c(option = "mako", direction = 1, na.value = "grey50", 
                         breaks = seq(0, -8, by = -2), labels = c("0", "2", "4", "6", ">8\nµg/m3"),
                         limits = c(-8, 1), oob = scales::squish) +
    labs(caption = "Difference in average indoor PM2.5 concentrations between no HVAC and HVAC scenarios",
         fill = NULL) +
    guides(fill = guide_colourbar(barwidth = 0.8, barheight = 9, label.position = "left")) +
    theme_void() + 
    coord_sf(clip = "off", xlim = c(-125, -114)) +
    theme(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(t = 10, b = 0)),
          plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
          legend.position = c(0.75, 0.80),
          legend.text = element_text(size = 7, colour = "grey20"),
          plot.title.position = "plot",
          panel.border = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

# # map indoor PM2.5 concentrations by tract for bay area
# p_inset <- ggplot(data = st_as_sf(df_tract), aes(fill = pm_25_diff_logic)) + 
#   geom_sf(colour = NA, size = 0.05, alpha = 0.95) +
#   scale_fill_viridis_c(option = "mako", direction = -1, na.value = "grey50") +
#   labs(title = NULL, subtitle = NULL, fill = NULL) +
#   guides(fill = "none") +
#   theme_void() + 
#   coord_sf(xlim = c(-122.6445, -121.5871), ylim =c(37.1897, 38.2033)) +
#   theme(panel.border = element_rect(colour = "grey20", fill = NA))
# 
# # built map with inset
# p1 <- p_map + 
#   annotate("rect", xmin = -122.6445, xmax = -121.5871, ymin = 37.1897, ymax = 38.2033, fill = NA, colour = "grey40", size = 0.2) +
#   #annotate("segment", x = c(-122.6445, -122.6445), xend = c(-123.88, -123.88), y = c(37.1897, 38.2033), yend = c(35.21, 38.305),
#   #         size = 0.20, colour = "grey40") +
#   inset_element(p_inset, left = -0.1, bottom = 0.3, right = 0.4, top = 0.6, clip = FALSE, on_top = FALSE) +
#   theme(plot.background = element_blank())


## TIMESERIES ##

# load simulation dataset
df_sim <- read_rds(here("data", "simulations", recent_df_hourly))

# filter to norcal counties
df_sim <- df_tract %>%
    select(geoid, county) %>%
    left_join(df_sim, ., by = "geoid") %>%
    filter(county %in% list_tracts$name)

# do timeseries plot of indoor PM2.5 concentrations
p2 <- df_sim %>%
    group_by(county,
             datetime = floor_date(datetime, "12 hours")) %>%
    summarise(pm_25_baseline = mean(pm_25_infiltration, na.rm = TRUE),
              pm_25_hvac = mean(pm_25_baseline_hvac, na.rm = TRUE),
              pm_25_logic = mean(pm_25_baseline_out, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with("pm_25"), names_to = "scenario", values_to = "pm_25") %>%
    mutate(scenario = str_remove_all(scenario, "pm_25_")) %>%
    ggplot(., aes(x = datetime, y = pm_25)) +
    geom_line(aes(colour = scenario, group = interaction(scenario, county)), alpha = 0.1, size = 0.4) +
    stat_summary(aes(group = scenario, colour = scenario), geom = "line", fun = mean, size = 1.0) +
    scale_x_datetime(expand = c(0, 0), date_breaks = "1 weeks", labels = label_date_short()) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(10, 40, by = 10), labels = c("10", "20", "30", "40\nµg/m3")) +
    scale_color_manual(values = c("baseline" = "#DA9781", "hvac" = "#779BBB", "logic" = "forestgreen"), 
                       labels = c("baseline" = "Infiltration", "hvac" = "baseline HVAC", "logic" = "Outdoor PM logic")) +
    labs(#caption = "Simulated indoor concentrations of PM2.5 for 3,721 census tracts in Northern California",
        x = NULL, y = NULL, colour = NULL) +
    guides(colour = guide_legend(label.position = "bottom", keywidth = 0.5, keyheight = 0.5, nrow = 1, label.hjust = 0.5)) +
    coord_cartesian(clip = "off", ylim = c(0, 50)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, colour = "grey20"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.5, 0.8),
          legend.direction = "horizontal")


# load outdoor pm2.5 dataset
df_epa <- read_rds(here("data", "US", "df_epa.rds"))

# get list of EPA stations
if(!file.exists(here("data", "US", "list_stations.rds"))) {
    list_stations <- airnow_loadAnnual(year = 2020, parameter = "PM2.5", dataDir = NULL) %>% # AIRNOW IS NOT WORKING IN R >= 4...
        monitor_subset(stateCodes = "CA") %>%
        monitor_extractMeta()
} else {
    list_stations <- read_rds(file = here("data", "US", "list_stations.rds"))
}

# reverse geocode stations to get county
if(!exists("list_stations_geo")) {
    list_stations_geo <- list_stations %>%
        tidygeocoder::reverse_geocode(lat = latitude, long = longitude, method = "osm", full_results = TRUE)
}

# drop stations in SoCal
list_stations_geo2 <- list_stations_geo %>%
    select("epa_station" = monitorID,
           county) %>%
    mutate(county = str_remove_all(str_to_lower(county), " county$")) %>%
    filter(county %in% list_tracts$name)

list_stations <- list_stations %>%
    rename(epa_station = monitorID)

# plot outdoor pm2.5 concentrations
p3 <- df_epa %>%
    filter(epa_station %in% list_stations$epa_station) %>%
    group_by(epa_station,
             datetime = floor_date(datetime, "6 hour")) %>%
    summarise(pm_25_out = mean(pm_25_out, na.rm = TRUE)) %>%
    ungroup() %>%
    left_join(., list_stations, by = "epa_station") %>%
    ggplot(., aes(x = datetime, y = pm_25_out)) +
    geom_line(aes(group = epa_station), colour = "grey80", alpha = 0.4, size = 0.2) +
    stat_summary(geom = "line", fun = mean, colour = "#94A5A0", size = 1.0) +
    scale_x_datetime(expand = c(0, 0), date_breaks = "1 weeks", labels = label_date_short()) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(100, 300, by = 100), labels = c("100", "200", "300\nµg/m3")) +
    labs(caption = "Outdoor concentrations of PM2.5 measured by 124 EPA stations in Northern California",
         x = NULL, y = NULL) +
    coord_cartesian(clip = "off", ylim = c(0, 330)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
          axis.text = element_text(size = 7, colour = "grey20"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.y.right = element_blank())


# build full plot
p1 + p2 + p3 +
    plot_layout(design = "AAABB
                        AAABB
                        AAACC") +
    plot_annotation(title = "PM2.5 Exposure from Wildfire Smoke",
                    subtitle = "PM2.5 concentrations for a synthetic dataset of California homes during Aug-Sep 2020 wildfire event") & 
    theme(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
          plot.background = element_rect(fill = "white", colour = "white"))

# save tract-level map
ggsave(here("products", "gfx", "plot_full_sim.png"), width = 3840, height = 1900, units = "px", dpi = 300)

# clean up
rm(p1, p2, p3, df_epa, list_stations, list_tracts)


#### Standard Error ####

# this is actually not making any sense as 
# Compute mean and standard deviation for each PM2.5 metric
tmp_df_summary <- df_sim %>%
    group_by(datetime) %>%
    summarise(
        mean_pm25_infiltration = mean(pm_25_infiltration, na.rm = TRUE),
        sd_pm25_infiltration = sd(pm_25_infiltration, na.rm = TRUE),
        mean_pm25_baseline_out = mean(pm_25_baseline_out, na.rm = TRUE),
        sd_pm25_baseline_out = sd(pm_25_baseline_out, na.rm = TRUE),
        mean_pm25_out = mean(pm_25_out, na.rm = TRUE),
        sd_pm25_out = sd(pm_25_out, na.rm = TRUE),
        .groups = "drop"
    )

# Reshape mean values to long format and rename variables
tmp_means <- tmp_df_summary %>%
    pivot_longer(
        cols = c(mean_pm25_infiltration, mean_pm25_baseline_out, mean_pm25_out),
        names_to = "variable",
        values_to = "mean"
    ) %>%
    mutate(variable = recode(variable,
                             "mean_pm25_infiltration" = "Sheltering Indoors",
                             "mean_pm25_baseline_out" = "Wildfire Mode",
                             "mean_pm25_out" = "Outdoor"
    ))

# Reshape standard deviations to long format and rename variables
tmp_sds <- tmp_df_summary %>%
    pivot_longer(
        cols = c(sd_pm25_infiltration, sd_pm25_baseline_out, sd_pm25_out),
        names_to = "variable",
        values_to = "sd"
    ) %>%
    mutate(variable = recode(variable,
                             "sd_pm25_infiltration" = "Sheltering Indoors",
                             "sd_pm25_baseline_out" = "Wildfire Mode",
                             "sd_pm25_out" = "Outdoor"
    ))

# Merge means and standard deviations
tmp_df_long <- left_join(tmp_means, tmp_sds, by = c("datetime", "variable"))

ggplot(tmp_df_long, aes(x = datetime, y = mean, color = variable)) +
    geom_line(size = 1) +  # Plot mean PM2.5 values
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = variable), 
                alpha = 0.2, color = NA) +  # Standard deviation shading
    labs(title = "PM2.5 Time Series with Standard Deviation",
         y = "PM2.5 Concentration (µg/m³)",
         x = "Date",
         color = "PM2.5 Type",
         fill = "PM2.5 Type") +
    theme_minimal()

ggsave(filename = here("products/gfx/fig8_stdev.png"), width = 3600, height = 1800, dpi = 300, units = "px")




#### FIGURE 8 ####

# Timeseries

ts_ts <- df_sim %>%
    group_by(county,
             datetime = floor_date(datetime, "1 hours")
    ) %>%
    summarise(pm_25_out = mean(pm_25_out, na.rm = TRUE),
              pm_25_baseline = mean(pm_25_infiltration, na.rm = TRUE),
              pm_25_hvac = mean(pm_25_baseline_hvac, na.rm = TRUE),
              pm_25_logic = mean(pm_25_baseline_out, na.rm = TRUE),
              pm_25_pac = mean(pm_25_pac, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(datetime = floor_date(datetime, "12 hours")) %>%
    summarise(pm_25_out = mean(pm_25_out, na.rm = TRUE),
              pm_25_baseline = mean(pm_25_baseline, na.rm = TRUE),
              pm_25_hvac = mean(pm_25_hvac, na.rm = TRUE),
              pm_25_logic = mean(pm_25_logic, na.rm = TRUE),
              pm_25_pac = mean(pm_25_pac, na.rm = TRUE)) %>%
    ggplot(., aes(x = datetime)) +
    geom_line(aes(y=pm_25_out), color = "black") +
    geom_line(aes(y=pm_25_baseline), color = "#FFB10A") +
    geom_line(aes(y=pm_25_pac), color = "#4AAEE8") +
    geom_line(aes(y=pm_25_hvac), color = "#5CC1A6") +
    geom_line(aes(y=pm_25_logic), color = "magenta") +
    #pivot_longer(cols = starts_with("pm_25"), names_to = "scenario", values_to = "pm_25") %>%
    #mutate(scenario = str_remove_all(scenario, "pm_25_")) %>%
    #ggplot(., aes(x = datetime, y = pm_25)) +
    #geom_line(aes(colour = scenario, group = interaction(scenario, county)), alpha = 0.1, size = 0.4) +
    #stat_summary(aes(group = scenario, colour = scenario), geom = "line", fun = mean, size = 1.0) +
    scale_x_datetime(expand = c(0, 0), date_breaks = "1 weeks", labels = label_date_short()) +
    #scale_y_continuous(expand = c(0, 0), breaks = seq(0, 120, by = 30), labels = c("", "30", "60", "90", "120\nµg/m3")) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 120, by = 30), labels = c("", "30", "60", "90", "120")) +
    #scale_color_manual(values = c("baseline" = "#DA9781", "hvac" = "#779BBB", "logic" = "forestgreen"), 
    #                   labels = c("baseline" = "Infiltration", "hvac" = "baseline HVAC", "logic" = "Outdoor PM logic")) +
    labs(#caption = "Simulated indoor concentrations of PM2.5 for 3,721 census tracts in Northern California",
        x = NULL, y = expression(PM[2.5] ~ "(µg/m³)"), #y = expression(atop(PM[2.5], "(µg/m³)")), 
        colour = NULL) +
    guides(colour = guide_legend(label.position = "bottom", keywidth = 0.5, keyheight = 0.5, nrow = 1, label.hjust = 0.5)) +
    coord_cartesian(clip = "off", ylim = c(0, 140)) +
    geom_segment(x = as.POSIXct("2020-09-13 12:00:00"), y = 115, xend = as.POSIXct("2020-09-13 12:00:00"), yend = 32,
                 arrow = arrow(length = unit(0.1, "cm")),
                 color = "#FFB10A", size = 0.4) +
    geom_segment(x = as.POSIXct("2020-09-13 12:00:00"), y = 28, xend = as.POSIXct("2020-09-13 12:00:00"), yend = 10,
                 arrow = arrow(length = unit(0.1, "cm")),
                 color = "magenta", size = 0.4) +
    geom_curve(aes(x = as.POSIXct("2020-09-02 12:00:00"), y = 92, 
                   xend = as.POSIXct("2020-09-12 10:00:00"), yend = 60),
               arrow = arrow(length = unit(0.1, "cm")), 
               curvature = 0.3,  # Adjust curvature: positive for upward, negative for downward
               color = "grey40", size = 0.4) +
    annotate("text", x = as.POSIXct("2020-09-02 00:00:00"), y = 98, 
             label = "Reduced Indoor\nConcentrations", hjust = 0.5, size = 3, color = "grey20") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
          plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, colour = "grey20"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.5, 0.8),
          legend.direction = "horizontal")

# Boxplots

df_sim %>%
    filter(datetime >= "2020-09-09", datetime <= "2020-09-16") %>%
    group_by(datetime = floor_date(datetime, "12 hours")) %>%
    summarise(pm_25_out = mean(pm_25_out, na.rm = TRUE),
              pm_25_infiltration = mean(pm_25_infiltration, na.rm = TRUE),
              pm_25_baseline_hvac = mean(pm_25_baseline_hvac, na.rm = TRUE),
              pm_25_baseline_out = mean(pm_25_baseline_out, na.rm = TRUE),
              pm_25_pac = mean(pm_25_pac, na.rm = TRUE)) %>%
    ungroup() %>%
    summarise(med_pm_25_out = median(pm_25_out, na.rm = TRUE),
              q3_pm_25_out = quantile(pm_25_out, 0.75, na.rm = TRUE),
              med_pm_25_infiltration = median(pm_25_infiltration, na.rm = TRUE),
              q3_pm_25_infiltration = quantile(pm_25_infiltration, 0.75, na.rm = TRUE),
              med_pm_25_pac = median(pm_25_pac, na.rm = TRUE),
              q3_pm_25_pac = quantile(pm_25_pac, 0.75, na.rm = TRUE),
              med_pm_25_baseline_hvac = median(pm_25_baseline_hvac, na.rm = TRUE),
              q3_pm_25_baseline_hvac = quantile(pm_25_baseline_hvac, 0.75, na.rm = TRUE),
              med_pm_25_baseline_out = median(pm_25_baseline_out, na.rm = TRUE),
              q3_pm_25_baseline_out = quantile(pm_25_baseline_out, 0.75, na.rm = TRUE)
    ) %>%
    glimpse()

ts_box <- df_sim %>%
    filter(datetime >= "2020-09-09", datetime <= "2020-09-16") %>%
    group_by(datetime = floor_date(datetime, "12 hours")) %>%
    summarise(pm_25_out = mean(pm_25_out, na.rm = TRUE),
              pm_25_infiltration = mean(pm_25_infiltration, na.rm = TRUE),
              pm_25_baseline_hvac = mean(pm_25_baseline_hvac, na.rm = TRUE),
              pm_25_baseline_out = mean(pm_25_baseline_out, na.rm = TRUE),
              pm_25_pac = mean(pm_25_pac, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(pm_25_out, pm_25_infiltration, pm_25_pac, pm_25_baseline_hvac, pm_25_baseline_out) %>%
    pivot_longer(cols = starts_with("pm"), names_to = "control_strategy", values_to = "concentration") %>%
    ggplot(., aes(x = as_factor(control_strategy) , y = concentration, color = control_strategy, fill = control_strategy)) +
    geom_boxplot(alpha = 0.6) +
    scale_color_manual(values = c("pm_25_out" = "black", 
                                  "pm_25_infiltration" = "#FFB10A", 
                                  "pm_25_baseline_hvac" = "#5CC1A6", 
                                  "pm_25_baseline_out" = "magenta", 
                                  "pm_25_pac" = "#4AAEE8")) +
    scale_fill_manual(values = c("pm_25_out" = "black", 
                                 "pm_25_infiltration" = "#FFB10A", 
                                 "pm_25_baseline_hvac" = "#5CC1A6", 
                                 "pm_25_baseline_out" = "magenta", 
                                 "pm_25_pac" = "#4AAEE8")) +
    guides(color = "none", fill = "none") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 120, by = 30), labels = c("", "30", "60", "90", "120\nµg/m3")) +
    coord_cartesian(clip = "off", ylim = c(0, 140)) +
    annotate("text", x = 6, y = 98.2, label = "Outdoors", hjust = 0, size = 3, color = "black") +
    annotate("text", x = 6, y = 24.5, label = "Sheltering Indoors", hjust = 0, size = 3, color = "#FFB10A") +
    annotate("text", x = 6, y = 20.1, label = "Portable air cleaner", hjust = 0, size = 3, color = "#4AAEE8") +
    annotate("text", x = 6, y = 14.7, label = "Air conditioning", hjust = 0, size = 3, color = "#5CC1A6") +
    annotate("text", x = 6, y = 7.26, label = "Wildfire mode", hjust = 0, size = 3, color = "magenta") +
    annotate("text", x = 1.2, y = 115, label = "101.2", hjust = 0, size = 3, color = "black") +
    annotate("text", x = 2.2, y = 32, label = "25.3", hjust = 0, size = 3, color = "#FFB10A") +
    annotate("text", x = 3.2, y = 28, label = "21", hjust = 0, size = 3, color = "#4AAEE8") +
    annotate("text", x = 4.2, y = 23, label = "14.4", hjust = 0, size = 3, color = "#5CC1A6") +
    annotate("text", x = 5.2, y = 12, label = "7.8", hjust = 0, size = 3, color = "magenta") +
    theme_minimal() +
    theme(
        plot.margin = margin(t = 5, r = 60, b = 5, l = 5),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )


design <- "
   11112
   11112
"

ts_ts + ts_box +
    plot_layout(design = design)

ggsave(filename = here("products","gfx","figure8.png"), 
       width = 3500, height = 1000, units = "px", dpi = 300)

#### ANIMATION ####

# reload simulation dataset
df_sim <- read_rds(here("data", "simulations", recent_df_hourly))

# load wildfire dataset and expand timeseries for counties
df_fires <- read_csv(here("data", "US", "df_fires.csv")) %>%
    mutate(across(where(is.character), str_to_lower)) %>%
    filter(end > start,
           acres > 300000) %>% # only include large fires
    mutate(start = force_tz(start, tzone = "America/Los_Angeles"),
           end = force_tz(end, tzone = "America/Los_Angeles"),
           datetime = map2(start, end, seq, by = 'day')) %>%
    unnest(cols = c(datetime)) %>%
    distinct(county, datetime) %>%
    group_by(county) %>%
    complete(datetime = seq.POSIXt(min(datetime), max(datetime), by = "hour")) %>%
    ungroup() %>%
    arrange(county, datetime) %>%
    mutate(fire = 1)

# select time of interest
df_map <- df_sim %>%
    mutate(datetime = with_tz(datetime, tzone = "America/Los_Angeles")) %>%
    filter(datetime >= "2020-09-07 12:00:00", 
           datetime <= "2020-09-10 12:00:00") %>%
    drop_na()

# aggregate to county level for speed
df_map <- df_tract %>%
    select(geoid, county) %>%
    left_join(df_map, ., by = "geoid") %>%
    mutate(pm_25_diff_logic = pm_25_baseline_out - pm_25_infiltration) %>%
    group_by(county,
             datetime) %>%
    summarise(pm_25_diff_logic = mean(pm_25_diff_logic, na.rm = TRUE)) %>%
    ungroup()

# add fires to dataframe
df_map <- df_map %>%
    left_join(., df_fires, by = c("county", "datetime")) %>%
    mutate(fire = replace_na(fire, 0))

# calculate running mean exposure
df_map <- df_map %>%
    group_by(county) %>%
    arrange(datetime, .by_group = TRUE) %>%
    mutate(samples = 1,
           pm_25_diff_run = cumsum(pm_25_diff_logic)/cumsum(samples)) %>% 
    ungroup() %>%
    select(-samples)

# add census geometry to dataframe
df_map <- df_county %>%
    select("county" = "NAME", 
           geometry) %>%
    mutate(county = str_to_lower(county)) %>%
    left_join(df_map, ., by = c("county"))

# define animated map
p_anim <- ggplot() + 
    geom_sf(data = st_as_sf(df_map), aes(fill = pm_25_diff_run, colour = fire > 0.5, size = fire > 0.5)) +
    #geom_sf(data = df_county, fill = NA, colour = "white", size = 0.03) +
    scale_colour_manual(values = c("white", "#E95635")) +
    scale_fill_gradient2(breaks = c(-50, -25, 0), labels = c("-50", "-25", "0 µg/m3"), midpoint = -12,
                         low = "#300A47", mid = "#3E6FA3", high = "#BCECF0", na.value = "grey50") +
    scale_size_manual(values = c(0.03, 0.12)) +
    labs(title = "Indoor Exposure to Wildfire Smoke",
         subtitle = "{format(frame_time, '%A %d %I%P\n%B')}",
         caption = "Difference in average indoor PM 2.5 concentrations between baseline and HVAC scenarios in Aug - Sep 2020\n for a synthetic dataset of houses with infiltration rates inferred from building characteristics",
         fill = NULL) +
    guides(fill = guide_colourbar(barwidth = 0.8, barheight = 9, label.position = "left"),
           colour = "none", size = "none") +
    transition_time(datetime) +
    theme_void() + 
    coord_sf(clip = "off", xlim = c(-125, -114)) +
    theme(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey20", face = "plain", hjust = 0.5, margin = margin(t = 12, b = -3)),
          plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
          legend.position = c(0.85, 0.8),
          legend.text = element_text(size = 7),
          plot.title.position = "plot",
          panel.border = element_blank(),
          plot.margin = margin(t = 0, r = 5, b = 0, l = 1, unit = "mm"))

# simplified map
# # Recreate the animation
# p_anim <- ggplot(df_map) +
#     geom_sf(aes(geometry = geometry, fill = pm_25_diff_logic)) +
#     transition_time(datetime) +
#     labs(title = "PM2.5 Difference Over Time: {frame_time}")
# 
# # Check the number of unique frames
# nframes <- length(unique(df_map$datetime))
# 
# # Generate animation
# map_anim <- animate(p_anim, 
#                     nframes = nframes, 
#                     fps = 12, 
#                     width = 6.5, height = 8.0, unit = "in", res = 150)
# 
# # Check the class
# class(map_anim)


# render animated map
map_anim <- animate(p_anim, nframes = length(unique(df_map$datetime)), 
                    fps = 12, start_pause = 24, end_pause = 48,
                    width = 6.5, height = 8.0, unit = "in", res = 150)

dim(df_map)  ## save animated map as gif
anim_save(here("products", "gfx", "anim_map_hourly.gif"), animation = map_anim)

# clean up
rm(p_anim, map_anim, df_fires, df_map)


## END ##