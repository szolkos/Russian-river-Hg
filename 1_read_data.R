#===========================================================================================================#
# 1_read_data.R
# Author: Scott Zolkos (sgzolkos@gmail.com) | January 14, 2022
# Background: Code for analyses and figures for Russian Arctic river long-term TSS and Hg
#===========================================================================================================#

  # BACKGROUND ####
  
    # Set working directory
      dir <- # SET WD HERE, e.g., "/Users/szolkos/Desktop/Russian_river_Hg_TSS_data/"
      setwd(dir)
      #dir <- "/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Data/for R/"
  
    # Load packages ####
      library(plotrix) # For plotting, labeling, axis and color scaling functions, etc.
      library(ggplot2) # For 'declaratively' creating graphics
      library(gdata) # To bind columns with different number of rows
      library(lubridate) # To convert from date to Julian day
      library(gridExtra)
      library(onewaytests) # for Welch's ANOVA for unequal variances
      library(ggpubr) # For multi-panel ggplot (e.g. 'ggarrange' function)
      library(fable) # For Dynamic regression model (DRM) (multivariate time-series analysis) | https://fable.tidyverts.org/index.html
      library(fpp3) # DRM examples from Forecasting Principles & Practice, 3rd ed. (Hyndman & Athanasopoulos 2021)
      library(feasts) # For DRM (Feature Extraction and Statistics for Time Series)
      library(tsfeatures) # For DRM (e.g., unitroot_kpss())
      library(magrittr) # For DRM (pipes)
      library(GGally) # For DRM (ggpairs)
      library(tibble) # For DRM
      library(plyr) # For DRM
      library(dplyr) # For DRM
      library(tidyr) # For DRM
      library(ggtext) # For PHg lit comparison
  
    # Packages used in development of code (listed here, but may or may not be needed)
      #library(DTK) # Dunnett-Tukey-Kramer Pairwise Multiple Comparison Test Adjusted for Unequal Variances and Unequal Sample Sizes
      #library(rstatix) # for repeated-measures ANOVA
      #library(nlme) # autocorrelation function
      #library(lme4) # For Linear Mixed Effects Regression
      #library(lmerTest) # For e.g. difflsmeans
      #library(car) # For Durbin-Watson Test, to evaluate correlation of residuals in a linear model
      #library(leaps) # regression subset selection
      #library(zyp) # For Kendall-Mann test of time series data, with Theil-Sen approach for resolving serial autocorrelation
      #library(Hmisc) # for rcorr function
      #library(corrplot) # for corrplot function
      #library(aTSA) # assess stationarity in time-series (e.g., adf.test())
      #library(foreign) # read .dbf file
  
    # Working in RStudio with GitHub:
      # https://cfss.uchicago.edu/setup/git-with-rstudio/
      # https://r-bio.github.io/intro-git-rstudio/
      
    # Annotated bibliography (references in alphabetical order)
      # Hirsch et al. (2010) - WRTDS (JAWRA, doi:/10.1111/j.1752-1688.2010.00482.x)
      # Hirsch et al. (2015) - WRTDS (USGS T&M 4-10: User Guide to Exploration and Graphics for RivEr Trends (EGRET) and dataRetrieval: R Packages for Hydrologic Data)
      # Hirsch et al. (2015) - WRTDS (EMS, doi:/10.1016/j.envsoft.2015.07.017)
      # Holmes et al. (2012) - for watershed areas and correcting spatiotemporal offsets (https://link.springer.com/article/10.1007/s12237-011-9386-6)
      # Hyndman & Athanasopoulos (2021) - for DRM (https://otexts.com/fpp3/dynamic.html | https://github.com/robjhyndman?tab=repositories)
      # Lehner et al. (2011) - GRanD reservoir database (https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/100125)
      # Pelletier et al. (2007) - (https://pubs.acs.org/doi/10.1021/acs.est.7b03929)
      # Zhang & Hirsch (2019) - WRTDS (WRR, doi:/10.1029/2019WR025338)
      # EGRET 3.0: http://usgs-r.github.io/EGRET/articles/Enhancements.html#generalized-flow-normalization-gfn-
    

  # FUNCTIONS #####
    
    # Plotting: no scientific format
      NotFancy <- function(l) {
      l <- format(l, scientific = FALSE)
      parse(text=l)
      }
    
    # Compile p values, R2 for linear model
      Lgnd <- function(mdl){
      if(class(mdl) != "lm") stop("Not an object of class 'lm' ")
      
      # Extract adjusted R2
      Rsq <- round(summary(mdl)$adj.r.squared,2)
      
      # Extract p-value for entire model
      f <- summary(mdl)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      if(p < 0.001){
        p <- 0.001
        lgnd <- c(bquote(italic(R)^2 == .(Rsq)), bquote(italic(p) < .(p)))
      }
      else if(0.001 < p & p < 0.01){
        p <- 0.01
        lgnd <- c(bquote(italic(R)^2 == .(Rsq)), bquote(italic(p) < .(p)))
      }
      else{
        p <- round(p,2)
        lgnd <- c(bquote(italic(R)^2 == .(Rsq)), bquote(italic(p) == .(p)))
      }
      return(lgnd)
    }
    
      
  # DATA PROCESSING & ANALYSIS ####
  
    # Definitions ####
      
      # Watershed areas (km2), total = 9,565,700 km2
        shed_area_on <- 55700 # Onega: 55700 # area from Arctic RIMS: http://rims.unh.edu/data/station/list.cgi?col=1
        shed_area_nd <- 348000 # NorthDvina: 348000 # area from Arctic RIMS: http://rims.unh.edu/data/station/list.cgi?col=1
        shed_area_me <- 564000 # Mezen: 564000 # area from Arctic RIMS: http://rims.unh.edu/data/station/list.cgi?col=1
        shed_area_pe <- 248000 # Pechora: 248000 # area from Arctic RIMS: http://rims.unh.edu/data/station/list.cgi?col=1
        shed_area_ob <- 2990000 # Ob: 2990000 # area from Holmes et al. 2012 E&C
        shed_area_ye <- 2400000 # Yenisey: 2400000 # area from Holmes et al. 2012 E&C
        shed_area_le <- 2430000 # Lena: 2430000 # area from Holmes et al. 2012 E&C
        shed_area_ko <- 530000 # Kolyma: 530000 # area from Holmes et al. 2012 E&C
        
      # Define years for which TSS and PHg were measured for each river
        on_tss_mmnts <- seq(1979,2018,1); on_phg_mmnts <- seq(1979,2002,1)
        nd_tss_mmnts <- seq(1979,2018,1); nd_phg_mmnts <- seq(1979,2002,1)
        me_tss_mmnts <- seq(1979,2018,1); me_phg_mmnts <- seq(1979,2002,1)
        pe_tss_mmnts <- seq(1979,2018,1); pe_phg_mmnts <- seq(1979,2002,1)
        ob_tss_mmnts <- seq(1974,2008,1); ob_phg_mmnts <- seq(1979,2008,1)
        ye_tss_mmnts <- seq(1974,2011,1); ye_phg_mmnts <- seq(1979,2011,1)
        le_tss_mmnts <- seq(2004,2018,1); le_phg_mmnts <- seq(2004,2011,1)
        ko_tss_mmnts <- seq(1980,1987,1); ko_phg_mmnts <- seq(1979,1987,1)
          
      # Set global parameters
        rus_rivers <- c("Onega","NorthDvina","Mezen","Pechora","Ob","Yenisey","Lena","Kolyma")
        rus_lgnd <- c("Onega","NorthDvina","Mezen","Pechora","Ob'","Yenisey","Lena","Kolyma")
        rus_lgnd2 <- c("Onega","NorthDvina","Mezen","Ob'","Yenisey","Lena","Kolyma")
        rus_riv_colors <- c("lightblue","orange","blue","gray","red","yellow","green1","purple3")
        rus_riv_colors2 <- c("lightblue","orange","blue","red","yellow","green1","purple3")
        rus_season_colors <- c("yellow","red","blue") # after Holmes et al. (2012)
        rus_season_lgnd <- c("Spring","Summer","Winter")
        rus_rivers_labels <- c("Onega","Northern\nDvina","Mezen","Pechora","Ob'","Yenisey","Lena","Kolyma")
      
      
    # Read data ####
      
      # Read data for TSS (mg/L) and Hg (ug/g)
        rus_hg <- read.csv(paste0(dir,"concentrations.csv"), header=T)
        
      # Discharge: m3/s; data obtained from Arctic-GRO website, A.I. Shiklomanov, and occasionally interpolated as needed following Methods detailed in this paper
        on_q <- read.csv(paste0(dir,"Q_Onega.csv"), header=T); on_q$Date <- as.Date(on_q$Date); on_q$Year <- year(on_q$Date); on_q <- droplevels(subset(on_q, select=c("River","Date","Year","Q")))
        nd_q <- read.csv(paste0(dir,"Q_NorthDvina.csv"), header=T); nd_q$Date <- as.Date(nd_q$Date); nd_q$Year <- year(nd_q$Date); nd_q <- droplevels(subset(nd_q, select=c("River","Date","Year","Q")))
        nd_q$River <- "NorthDvina"
        me_q <- read.csv(paste0(dir,"Q_Mezen.csv"), header=T); me_q$Date <- as.Date(me_q$Date); me_q$Year <- year(me_q$Date); me_q <- droplevels(subset(me_q, select=c("River","Date","Year","Q")))
        pe_q <- read.csv(paste0(dir,"Q_Pechora.csv"), header=T); pe_q$Date <- as.Date(pe_q$Date); pe_q$Year <- year(pe_q$Date); pe_q <- droplevels(subset(pe_q, select=c("River","Date","Year","Q")))
        ob_q <- read.csv(paste0(dir,"Q_Ob.csv"), header=T); ob_q$Date <- as.Date(ob_q$Date); ob_q$Year <- year(ob_q$Date); ob_q <- droplevels(subset(ob_q, select=c("River","Date","Year","Q")))
        ye_q <- read.csv(paste0(dir,"Q_Yenisey.csv"), header=T); ye_q$Date <- as.Date(ye_q$Date); ye_q$Year <- year(ye_q$Date); ye_q <- droplevels(subset(ye_q, select=c("River","Date","Year","Q")))
        le_q <- read.csv(paste0(dir,"Q_Lena.csv"), header=T); le_q$Date <- as.Date(le_q$Date); le_q$Year <- year(le_q$Date); le_q <- droplevels(subset(le_q, select=c("River","Date","Year","Q")))
        ko_q <- read.csv(paste0(dir,"Q_Kolyma.csv"), header=T); ko_q$Date <- as.Date(ko_q$Date); ko_q$Year <- year(ko_q$Date); ko_q <- droplevels(subset(ko_q, select=c("River","Date","Year","Q")))
        