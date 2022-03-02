#===========================================================================================================#
# 6b_climate.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | December 27, 2021
# Background: Trends in mean annual air temperature (MAAT) and total annual precipitation (TAP) by watershed
#===========================================================================================================#

  # BACKGROUND ####

    # Parameters & Units ####
      # MAAT = ˚C | Mean annual air temperature
      # TAP = m | Total annual precipitation
          

  # FUNCTIONS ####
  
    era_clim <- function(River, TempOrPrecip){

      # Subset data
        if(TempOrPrecip=="Temp"){
          era_sub <- droplevels(subset(era, era$river==River, select=c("year", "river", "maat", "maat_sd")))
          ylab <- "Mean annual air temperature (˚C)"
          lms <- c(-15,5)
          brks <- seq(-15,5,5)
        }
        if(TempOrPrecip=="Precip"){
          era_sub <- droplevels(subset(era, era$river==River, select=c("year", "river", "tap", "tap_sd")))
          era_sub$tap <- era_sub$tap*1000
          era_sub$tap_sd <- era_sub$tap_sd*1000
          ylab <- "Total annual precipitation (mm)"
          lms <- c(0,1300)
          brks <- seq(0,1200,200)
        }
        names(era_sub) <- c("year", "river", "mean", "sd")
        
      # Plotting conditionals
        if(River=="Onega"){clr <- "lightblue"}
        if(River=="NorthDvina"){clr <- "orange"}
        if(River=="Mezen"){clr <- "blue"}
        if(River=="Pechora"){clr <- "gray"}
        if(River=="Ob"){clr <- "red"}
        if(River=="Yenisey"){clr <- "yellow"}
        if(River=="Lena"){clr <- "green1"}
        if(River=="Kolyma"){clr <- "purple3"}
        
      # Plot
        ggplot(data=era_sub, aes(y=mean, x=year)) +
          # Scales
            scale_y_continuous(limits=lms, breaks=brks) +
          # Themes
            theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.background=element_rect(fill="white"),
                  panel.border=element_rect(colour="black", fill=NA, size=1),
                  plot.background=element_rect(colour="white", size=1),
                  axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
                  axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0)),
                  #axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10)),
                  axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
                  axis.title.x=element_blank(),
                  #axis.title.x=element_text(margin=margin(t=-10, r=0, b=0, l=0)),
                  text=element_text(size=19),
                  legend.position="none") +
          # Data
            stat_summary(data=era_sub, aes(ymin=mean-sd, ymax=mean+sd), fill=clr, fun.data=mean, geom="ribbon", alpha=0.5) +
            geom_path(color="black", size=1) +
            geom_point(data=era_sub, aes(y=mean, x=year), pch=21, fill=clr, col="black", size=4, stroke=1, alpha=1) +
          # Axis labels
            labs(title=River, y=ylab)
      }
    
    climate_trends <- function(TSSorPHg){

      # Store durations of sampling
        if(TSSorPHg == "TSS"){
          on_len <- c(1980,1992)
          nd_len <- c(1979,2001)
          me_len <- c(1980,2018)
          pe_len <- c(1983,2018)
          ob_len <- c(1979,2008)
          ye_len <- c(1980,2011)
          le_len <- c(2004,2018)
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
      
      # Create dataframe
        Names <- c("river", "param", "t1", "t2", "val_t1", "val_t2", "change_abs", "change_pcnt", "MK_sig", "MK_ann_trend")
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
          # Data
            era_sub <- droplevels(subset(era, era$river==i & era$year >= len[1] & era$year <= len[2]))
            era_sub$tap <- era_sub$tap*1000 # to convert m to mm for precip
          # Populate df
            df[j,1] <- i
            df[j,2] <- "MAAT"
            df[j,3] <- start <- len[1] # start year
            df[j,4] <- end <- len[2] # end year
            df[j,5] <- val_t1 <- round(subset(era_sub, era_sub$year == start)$maat,2) # MAAT (˚C) at t1
            df[j,6] <- val_t2 <- round(subset(era_sub, era_sub$year == end)$maat,2) # MAAT at t2
            df[j,7] <- val_t2-val_t1 # absolute change
            df[j,8] <- round(((val_t2-val_t1)/abs(val_t1))*100,1) # relative (%) change
            df[j,9] <- round(zyp.trend.vector(era_sub$maat)[6],3) # Mann-Kendall significance
            df[j,10] <- round(zyp.trend.vector(era_sub$maat)[2],3) # trend (annual change)
            df[j+8,1] <- i
            df[j+8,2] <- "TAP"
            df[j+8,3] <- start <- len[1] # start year
            df[j+8,4] <- end <- len[2] # end year
            df[j+8,5] <- val_t1 <- round(subset(era_sub, era_sub$year == start)$tap,1) # TAP (mm) at t1
            df[j+8,6] <- val_t2 <- round(subset(era_sub, era_sub$year == end)$tap,1) # TAP at t2
            df[j+8,7] <- val_t2-val_t1 # absolute change
            df[j+8,8] <- round(((val_t2-val_t1)/abs(val_t1))*100,1) # relative (%) change
            df[j+8,9] <- round(zyp.trend.vector(era_sub$tap)[6],3) # Mann-Kendall significance
            df[j+8,10] <- round(zyp.trend.vector(era_sub$tap)[2],3) # trend (annual change)
          # Iterative
            j <- j+1
      }
      
      return(df)
    }

      
  # DATA PROCESSING AND ANALYSIS ####
    
    # Data ####
  
      # ERA5 climate data: mean annual air temperature (MAAT) and total annual precipitation (TAP), 1979-2000 from Google Earth Engine
        era <- read.csv(paste0(dir,"era5_temp_precip_1979_2020.csv"), header=T)
        names(era) <- c("year","river","maat","maat_sd","tap","tap_sd")
        era$river <- factor(era$river, levels=rus_rivers)
    
    # Statistics ####
        
      # Changes in MAAT and TAP by watershed
        clim_trends_df <- climate_trends("PHg")
      # Optionally export
        #write.csv(clim_trends_df, paste0(dir,"clim_trends_df.csv"), row.names=F)
          
    # Figures ####
    
      # MAAT or TAP by year w/ SD as shaded region- Selected river (export as 5" x 6" landscape PDF)
        era_clim("Onega", "Temp")
        for(i in rus_rivers){plot(era_clim(i, "Temp"))}
        
          
  # Plotting ####
          
    # Example: MAAT by year- All rivers
  
        ggplot(data=era, aes(y=maat, x=year, fill=river)) +
          # Scales
            scale_y_continuous(limits=c(-15,5), breaks=seq(-15,5,5)) +
            # Themes
            theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank(),
                  plot.background=element_rect(colour="white", size=1),
                  panel.border=element_rect(colour="black", fill=NA, size=1),
                  text=element_text(size=19)) +
            theme(axis.text.y=element_text(angle=0, hjust=0.5, colour="black")) +
            theme(axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0))) +
            theme(axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10))) +
            theme(axis.text.x=element_text(angle=0, hjust=0.5, colour="black")) +
            theme(axis.title.x=element_blank()) +
            #theme(axis.title.x=element_text(margin=margin(t=-10, r=0, b=0, l=0))) +
            theme(plot.background=element_rect(fill='white')) +
            theme(panel.background=element_rect(fill='white')) +
            theme(legend.position="none") +
          # Data
            #stat_summary(geom="ribbon", aes(ymin=maat-maat_sd, ymax=maat+maat_sd), fill=rus_riv_colors[era$river], alpha=0.5) +
            #stat_summary(data=era, geom="ribbon", fun=mean, fun.min=era$maat-era$maat_sd, fun.max=era$maat+era$maat_sd, fill="gray40", alpha=0.5) + #rus_riv_colors[era$river]
            #stat_summary(data=era, aes(fill=river, ymin=maat-maat_sd, ymax=maat+maat_sd), fun.data=mean, geom="ribbon", alpha=0.5) +
            stat_summary(data=era, aes(ymin=maat-maat_sd, ymax=maat+maat_sd), fun.data=mean, geom="ribbon", alpha=0.5) +
            geom_path(color=rus_riv_colors[era$river], size=1) +
            #geom_errorbar(data=era, aes(y=maat, x=year, ymin=maat-maat_sd, ymax=maat+maat_sd), color=rus_riv_colors[era$river], height=1, position=position_dodge()) +
            geom_point(data=era, aes(y=maat, x=year), pch=21, fill=rus_riv_colors[era$river], col="black", size=4, stroke=1, alpha=1) +
          # Axis labels
            labs(y="Mean annual air temperature (˚C)", x="")
          
  