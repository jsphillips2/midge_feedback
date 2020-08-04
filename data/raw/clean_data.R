#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(readxl)
library(lubridate)

# set import path
import <- "data/raw/jp_midgeprod_main_20sep17.xlsx"

# examine sheets
sheets <- excel_sheets(import)

# import data
cores <- read_excel(import, sheet="Cores")
inc <- read_excel(import, sheet="Incubations")
hobo <- read_excel(import, sheet="HOBO")
adults = read_excel(import, sheet="Adults Check")

# convert to dataes dates
inc <- inc %>% mutate(sampledate = as.Date(sampledate))
hobo <- hobo %>% mutate(sampledate = as.Date(sampledate))






#==========
#========== HOBO: Light and Temperature
#==========

# export csv
# write_csv(hobo, "data/clean/hobo.csv")

# extract incubation interval
inc_int <- cores %>%
  select(-comments) %>%
  left_join(inc) %>%
  group_by(rack, sampledate, shade) %>%
  summarize(
    start = max(time_start),
    end = min(time_end)
  ) 

# expand hobo data by light and dark (necessary for mathcing incubations)
hobo_shade <- hobo %>%
  tidyr::expand(nesting(rack, sampledate, sampletime, temperature, light_intensity),
         shade = c("light","dark")) 

# join hobo_shade and inc_int, selecting data with matching dates 
# define incubation interval
hobo_inc <- hobo_shade %>%
  inner_join(inc_int) %>%
  mutate(
    start_end = interval(start, end)
  )

# check hobo_inc dimensions
# racks 1-4 should have 48 measurements for each of light dark on each day (96)
# racks 5,6 should have 96 measurements for each of light dark on each day (192)
hobo_inc %>%
  group_by(sampledate, rack) %>%
  summarize(samples = length(sampletime)) 

# check hobo_inc values
# make sure dates and times make sense
hobo_inc %>% summary
hobo_inc %>% 
  select(rack, shade, sampledate, start_end) %>% 
  as.data.frame %>% head(100)

# add logical to indicate which hobo readings are within the incubation period
hobo_inc2 <- hobo_inc %>%
  mutate(
    within_inc = sampletime %within% start_end
  )

# check values (dplyr doesn't like Period and Interval classes...)
hobo_inc2[hobo_inc2$rack==1&hobo_inc2$within_inc==F,
          c("sampletime","start_end","within_inc")] %>% head(100)
hobo_inc2[hobo_inc2$rack==1&hobo_inc2$within_inc==T,
          c("sampletime","start_end","within_inc")] %>% head(100)


# select hobo values within incubation window and aggregate
# convert light intensity to PAR (divide by 54)
hobo_summary <- hobo_inc2[hobo_inc2$within_inc==T,
                         c("rack","shade","sampledate","temperature","light_intensity")] %>%
  group_by(rack, shade, sampledate) %>%
  summarize(
    temperature = mean(temperature),
    par = mean(light_intensity/54)
  )

hobo_summary %>% as.data.frame %>% head(24)





#==========
#========== Incubations
#==========

# calculate change in do concentration
# divide mg/L by 1000 to convert to cubic cm
# multiply by water column depth to yield mass of oxygen per area
# multiply mg/cm^2 by 10000 to convert to mg/m^2
# divide by 1000 to convert to g/m^2
# divide by time interval
# divide midge_treat by cross_section_area to convert to #/cm^2
# multiply midge_treat by 10000 to convert to #/m^2

inc_proc <- cores %>%
  select(-comments) %>%
  left_join(inc) %>%
  mutate(o2_flux = (do_end - do_start)*column_depth/(100*duration),
    midge_density = round(10000*midge_treat/cross_section_area)) %>%
  select(core, rack, trial, sampledate, shade, site, midge_density, midge_treat, o2_flux)

# calculate NEP, ER, and GPP
met_wide <- inc_proc %>%
  spread(shade, o2_flux) %>%
  rename(nep = light, er = dark) %>%
  filter(er < 0) %>%
  mutate(gpp = nep - er)

# convert to long form and join with HOBO data
# convert trial to binary indicator
met_long <- met_wide %>%
  gather(type, met, gpp, er, nep) %>%
  mutate(shade = ifelse(type=="er","dark","light"),
         trial = trial - 1) %>%
  full_join(hobo_summary)

# export
# met_long %>% write_csv("data/clean/met_long.csv")





#==========
#========== Adults
#==========

# combine core and incubation information
adults = 
  cores %>%
  select(-comments) %>%
  left_join(adults) %>%
  # clean up variables
  # divide midge_treat by cross_section_area to convert to #/cm^2
  # multiply midge_treat by 10000 to convert to #/m^2
  mutate(
    sampledate = as.Date(sampledate),
    site = factor(site),
    midge_density = round(10000*midge_treat/cross_section_area)
  ) %>%
  # group by core sort
  group_by(core) %>%
  arrange(core, sampledate) %>%
  mutate(
    # calculate cummulative emergence on a given date
    cumem = cumsum(adults),
    
    # calculate porportion of midges emerged on a given date
    prop = ifelse(midge_treat > 0, cumem/midge_treat,NA),
    
    # calculate number of days in sampling period
    days = c(1, diff(sampledate)),
    
    # calculate rate at which adults emerged
    adults_rate = adults/days
  ) %>%
  ungroup()

# export
# adults %>% write_csv("data/clean/adults.csv")







