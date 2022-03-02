#===========================================================================================================#
# 4_Hg_dep.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | January 18, 2022
# Background: Read in data for Hg deposition from GEOS-Chem.
#===========================================================================================================#

  # BACKGROUND ####

    # Parameters & Units ####
  
      # regional_n = Hg deposition (Mg/y) from regional (EU + USSR) sources, lagged by n years
      # global_n = Hg deposition (Mg/y) from global sources, lagged by n years

  
  # FUNCTIONS ####
  
    Hg_dep_trends <- function(TSSorPHg){
      
      # Store durations
        if(TSSorPHg == "TSS"){
          on_len <- c(1980,1992)
          nd_len <- c(1979,2001)
          me_len <- c(1980,2015)
          pe_len <- c(1983,2015)
          ob_len <- c(1974,2008)
          ye_len <- c(1980,2011)
          le_len <- c(2004,2015)
          ko_len <- c(1982,1987)
        }
        if(TSSorPHg == "PHg"){
          on_len <- c(1980,1992)
          nd_len <- c(1980,2001)
          me_len <- c(1980,2002)
          pe_len <- c(1983,2002)
          ob_len <- c(1979,2008)
          ye_len <- c(1980,2011)
          le_len <- c(2004,2011)
          ko_len <- c(1982,1987)
        }
      
      # Data for trends analysis of contemporary Hg deposition
        Hg_dep_trend <- droplevels(subset(Hg_dep, select=c("river", "year", "regional_0", "global_0")))
      
      # Create dataframe
        Names <- c("river", "t1", "t2", "val_t1", "val_t2", "change_abs", "change_pcnt")
        df <- as.data.frame(matrix(ncol=length(Names), nrow=length(rus_rivers)*2))
        names(df) <- Names
        j=1
      # Compile Hg deposition trends
        for(i in rus_rivers){
          # Conditional for duration, based on river
            if(i=="Onega"){len <- on_len}
            if(i=="NorthDvina"){len <- nd_len}
            if(i=="Mezen"){len <- me_len}
            if(i=="Pechora"){len <- pe_len}
            if(i=="Ob"){len <- ob_len}
            if(i=="Yenisey"){len <- ye_len}
            if(i=="Lena"){len <- le_len}
            if(i=="Kolyma"){len <- ko_len}
          # Populate df
            df[j,1] <- i
            df[j,2] <- start <- len[1] # start year
            df[j,3] <- end <- len[2] # end year
            df[j,4] <- val_t1 <- subset(Hg_dep_trend, Hg_dep_trend$river == i & Hg_dep_trend$year == start)$global_0 # contemporary global dep at t1
            df[j,5] <- val_t2 <- subset(Hg_dep_trend, Hg_dep_trend$river == i & Hg_dep_trend$year == end)$global_0 # contemporary global dep at t2
            df[j,6] <- val_t2-val_t1 # absolute change
            df[j,7] <- ((val_t2-val_t1)/abs(val_t1))*100 # relative (%) change
            df[j+8,1] <- i
            df[j+8,2] <- start <- len[1] # start year
            df[j+8,3] <- end <- len[2] # end year
            df[j+8,4] <- val_t1 <- subset(Hg_dep_trend, Hg_dep_trend$river == i & Hg_dep_trend$year == start)$regional_0 # contemporary regional dep at t1
            df[j+8,5] <- val_t2 <- subset(Hg_dep_trend, Hg_dep_trend$river == i & Hg_dep_trend$year == end)$regional_0 # contemporary regional dep at t2
            df[j+8,6] <- val_t2-val_t1 # absolute change
            df[j+8,7] <- ((val_t2-val_t1)/abs(val_t1))*100 # relative (%) change
          # Iterative
            j <- j+1
        }
      
      return(df)
    }

    
  # DATA PROCESSING AND ANALYSIS ####
    
    # Data ####
      
      # Read data
        dir <- "/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Data/for R/"
        Hg_dep <- read.csv(paste0(dir,"GEOS_Chem_Hg_dep.csv"), header=T)
        
      # Compile Hg deposition trends
        phg_dep_trends <- Hg_dep_trends("PHg")
        phg_dep_trends
        
      # Optionally export
        dir <- "/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Data/from R/"
        #write.csv(phg_dep_trends, paste0(dir,"phg_dep_trends.csv"), row.names=F)
        