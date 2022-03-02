#===========================================================================================================#
# 5_reservoirs.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | December 27, 2021
# Background: Data from GRanD on reservoirs in Russian river watersheds
#===========================================================================================================#

  # DATA PROCESSING AND ANALYSIS ####

    # Data ####

      # Read reservoir data (from GRanD (Lehnher et al., 2011)
        grand_res_Russia <- read.csv(paste0(dir,"GRanD_Russia.csv"), header=T)
        
      # Explore data: subset by river
        river_reservoir <- "Vilyui"
        droplevels(subset(grand_res, grand_res$RIVER==river_reservoir, select=c("GRAND_ID","LAT_DD","LONG_DD","YEAR","RIVER","DAM_NAME","RES_NAME","AREA_SKM","AREA_REP","CAP_REP"))) 
        
      # Read and compile data on reservoir capacities (km3)
        res_cum_cap <- read.csv(paste0(dir,"res_cum_cap_w_lags.csv"), header=T) # Capacity of each reservoir and cumulative capacities for each river through time, with 5 year lags from 0-25 years
        res_cum_cap <- droplevels(subset(res_cum_cap, select=c("year", "river", "res_cum_cap_20", "res_cum_cap_25"), res_cum_cap$year >= 1979))
        