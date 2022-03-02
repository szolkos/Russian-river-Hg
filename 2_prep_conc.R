#===========================================================================================================#
# 2_prep_conc.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | December 23, 2021
# Background: Prepare concentration data:
#             (1) Calculate mean of replicates;
#             (2) Correct for C-Q spatiotemporal offsets;
#             (3) Calculate PHg.
#===========================================================================================================#

  # BACKGROUND ####

    # Parameters & units ####
      # Raw data (in 'conc_df')
        # TSS = mg/L | concentration of total suspended solids (TSS)
        # Hg = µg/g | concentration of mercury in TSS
        # PHg = ng/L | Hg*TSS
  

  # FUNCTIONS ####
  
    # Subset and derive annual [PHg] and [TSS] ± SE for desired river
      ann_conc <- function(river){
        # Subset data
          cp_df <- droplevels(subset(conc_df, select=c("Date","Hg","TSS","Season"), conc_df$River==river))
          cp_df$ratio <- cp_df$Hg/cp_df$TSS
          cp_df$Year <- year(cp_df$Date)
          cp_hg_mean <- as.data.frame(as.matrix(tapply(cp_df$Hg, cp_df$Year, FUN=mean, na.rm=T)))
          cp_hg_se <- as.data.frame(as.matrix(tapply(cp_df$Hg, cp_df$Year, FUN=std.error, na.rm=T)))
          cp_tss_mean <- as.data.frame(as.matrix(tapply(cp_df$TSS, cp_df$Year, FUN=mean, na.rm=T)))
          cp_tss_se <- as.data.frame(as.matrix(tapply(cp_df$TSS, cp_df$Year, FUN=std.error, na.rm=T)))
          cp_ratio_mean <- as.data.frame(as.matrix(tapply(cp_df$ratio, cp_df$Year, FUN=mean, na.rm=T)))
          cp_ratio_se <- as.data.frame(as.matrix(tapply(cp_df$ratio, cp_df$Year, FUN=std.error, na.rm=T)))
        # Bind
          cp_hg_tss <- cbind(river,row.names(cp_hg_mean),cp_hg_mean,cp_hg_se,cp_tss_mean,cp_tss_se,cp_ratio_mean,cp_ratio_se)
          names(cp_hg_tss) <- c("River","Year","Hg_mean","Hg_SE","TSS_mean","TSS_SE","Ratio_mean","Ratio_SE")
          cp_hg_tss$Hg_mean[cp_hg_tss$Hg_mean==0.00] <- NA
          cp_hg_tss$TSS_mean[cp_hg_tss$TSS_mean==0.00] <- NA
          cp_hg_tss$Ratio_mean[cp_hg_tss$Ratio_mean==0.00] <- NA
          cp_hg_tss$Year <- as.numeric(levels(cp_hg_tss$Year))
        # Replace NaN w/ NA
          cp_hg_tss$Hg_mean[cp_hg_tss$Hg_mean=="NaN"] <- NA
          cp_hg_tss$Ratio_mean[cp_hg_tss$Ratio_mean=="NaN"] <- NA
        # Return df
        return(cp_hg_tss)
      }


  # DATA PROCESSING AND ANALYSIS ####

    # 1. Calculate mean TSS and PHg for days with replicate samples ####
    
      # Explore: For each river, how many days with replicates?
        # TSS
          for(i in rus_rivers){
            pre_tss_agg <- droplevels(subset(rus_hg, rus_hg$River==i))
            tss_agg <- aggregate(pre_tss_agg$TSS, list(Date=pre_tss_agg$Date, River=pre_tss_agg$River), data=rus_hg, FUN=mean)
            head(pre_tss_agg); head(tss_agg)
            dim(pre_tss_agg); dim(tss_agg)
            print(paste0(i," (TSS): total n = ", nrow(pre_tss_agg), ", duplicate n = ", nrow(pre_tss_agg)-nrow(tss_agg), ",  net n = ", nrow(tss_agg)))
          }
        # PHg
          for(i in rus_rivers){
            pre_hg_agg <- droplevels(subset(rus_hg, rus_hg$River==i))
            hg_agg <- aggregate(pre_hg_agg$Hg, list(Date=pre_hg_agg$Date, River=pre_hg_agg$River), data=rus_hg, FUN=mean)
            head(pre_hg_agg); head(hg_agg)
            dim(pre_hg_agg); dim(hg_agg)
            print(paste0(i," (Hg): total n = ", nrow(pre_hg_agg), ", duplicate n = ", nrow(pre_hg_agg)-nrow(hg_agg), ",  net n = ", nrow(hg_agg)))
          }
      
      # Calculate means of replicates
        # TSS
          tss_agg <- aggregate(rus_hg$TSS, list(Date=rus_hg$Date, River=rus_hg$River), data=rus_hg, FUN=mean)
          names(tss_agg) <- c("Date","River","TSS"); head(tss_agg); dim(tss_agg)
        # PHg
          hg_tss_agg <- aggregate(rus_hg$Hg, list(Date=rus_hg$Date, River=rus_hg$River), data=rus_hg, FUN=mean)
          names(hg_tss_agg) <- c("Date","River","Hg"); head(hg_tss_agg); dim(hg_tss_agg)
        
      # Merge TSS and Hg by year and River
        conc_df <- merge(x=tss_agg, y=hg_tss_agg, by=c("River","Date"), all.x=T, all.y=T)
        head(conc_df)
  
        
    # 2. Correct constituent sampling date for spatiotemporal offset with discharge ####
  
      # Do this BEFORE merging C and Q data
      # Following Holmes et al. (2012), add offset to concentration mmnt date, to match with Q measurement date
      # Distance (km) which constituent measurements downstream for, assuming flow of 1.5 m/s: (e.g. Lena = (520 km * 1000 m/km) ÷ (1.5 m/s * 86400 s/d) =  d)
      # Rule: if C downstream of Q, substract offset from C sampling day so it becomes same as Q sample date
      # Offsets:
        # Mezen 100 km (chem collected downstream of Q [Q mmnt at Malonisogorskoye]) = (100*1000)/(1.5*86400) = 0.77 => -1 day
        # Pechora 315 km (chem collected downstream of Q [Q mmnt at Ust-Tsilma]) = (315*1000)/(1.5*86400) = 2.43 => -2 days
        # Yenisey 250 km (chem collected downstream of Q [Q mmnt at Igarka]) = (250*1000)/(1.5*86400) = 1.93 => -2 days
        # Lena 520 km (chem collected upstream of Q [Q mmnt at Kyusyur]) = (520*1000)/(1.5*86400) = 4.01 => +4 days
        # Kolyma 160 km (chem collected downstream of Q [Q mmnt at Kolymskoye]) = (160*1000)/(1.5*86400) = 1.23 => -1 day
        
      # Store dates as date
        conc_df$Date <- as.Date(conc_df$Date)
      # Apply corrections
        conc_df$Date[conc_df$River=="Mezen"] <- conc_df$Date[conc_df$River=="Mezen"]-1
        conc_df$Date[conc_df$River=="Pechora"] <- conc_df$Date[conc_df$River=="Pechora"]-2
        conc_df$Date[conc_df$River=="Yenisey"] <- conc_df$Date[conc_df$River=="Yenisey"]-2
        conc_df$Date[conc_df$River=="Lena"] <- conc_df$Date[conc_df$River=="Lena"]+4
        conc_df$Date[conc_df$River=="Kolyma"] <- conc_df$Date[conc_df$River=="Kolyma"]-1
        head(conc_df); dim(conc_df); length(na.omit(conc_df$Hg)); length(na.omit(conc_df$TSS))
  
        
    # 3. QA/QC ####
      # Remove Hg values below DL
        conc_df$HgperTSS <- conc_df$Hg/conc_df$TSS
        conc_df$Hg[conc_df$HgperTSS <= 0.005] <- NA
        conc_df$HgperTSS <- NULL
        head(conc_df); dim(conc_df); length(na.omit(conc_df$Hg)); length(na.omit(conc_df$TSS))
      # Omit anomalous Mezen Hg value (75 ng/L) on 1980-08-18
        conc_df$Hg[conc_df$Hg>50] <- NA
        head(conc_df); dim(conc_df); length(na.omit(conc_df$Hg)); length(na.omit(conc_df$TSS))
      # Remove TSS values below DL
        conc_df$TSS[conc_df$TSS <= 3] <- NA
        head(conc_df); dim(conc_df); length(na.omit(conc_df$Hg)); length(na.omit(conc_df$TSS))
      # Omit values from 1979 for Onega, NorthDvina, Mezen, where n=1 and in December, thus skewing results
        conc_df$TSS[conc_df$River=="Onega" & year(conc_df$Date)==1979] <- NA
        conc_df$Hg[conc_df$River=="Onega" & year(conc_df$Date)==1979] <- NA
        conc_df$Hg[conc_df$River=="NorthDvina" & year(conc_df$Date)==1979] <- NA
        conc_df$TSS[conc_df$River=="Mezen" & year(conc_df$Date)==1979] <- NA
        conc_df$Hg[conc_df$River=="Mezen" & year(conc_df$Date)==1979] <- NA
      # Calculate concentration of particulate Hg (PHg; i.e., Hg in TSS)
        conc_df$pHg <- conc_df$Hg*conc_df$TSS
        head(conc_df); dim(conc_df); length(na.omit(conc_df$Hg)); length(na.omit(conc_df$TSS)); length(na.omit(conc_df$pHg))
      # Add seasons
        conc_df$Season <- as.numeric(substr(conc_df$Date,6,7))
        conc_df$Season[conc_df$Season==5 | conc_df$Season==6] <- "Spring"
        conc_df$Season[conc_df$Season==7 | conc_df$Season==8 | conc_df$Season==9 | conc_df$Season==10] <- "Summer"
        conc_df$Season[conc_df$Season==11 | conc_df$Season==12 | conc_df$Season==1 | conc_df$Season==2 | conc_df$Season==3 | conc_df$Season==4] <- "Winter"
        conc_df$Season <- factor(conc_df$Season, levels=c("Spring","Summer","Winter"))
        conc_df$River <- factor(conc_df$River, levels=rus_rivers)
        
  
    # 4. Calculations ####
        
        dim(na.omit(subset(rus_hg, rus_hg$River=="Yenisey", select=c("Date","River","TSS","Hg"))))
        dim(na.omit(subset(conc_df, conc_df$River=="Yenisey", select=c("Date","River","pHg"))))
        
      # Annual [PHg] and [TSS] ± SE
        # For desired river
          ann_conc("Yenisey")
        # All rivers
          ann_conc_df <- rbind(ann_conc("Onega"), ann_conc("NorthDvina"), ann_conc("Mezen"), ann_conc("Pechora"), ann_conc("Ob"), ann_conc("Yenisey"), ann_conc("Lena"), ann_conc("Kolyma"))
        
      # Summaries
        # Select river and parameter
          riv <- "Ob" # Onega, NorthDvina, Mezen, Pechora, Ob, Yenisey, Lena, Kolyma
          param <- "pHg" # TSS, Hg, pHg
        # How many measurements?
          length(na.omit(subset(conc_df, conc_df$River==riv, select=c("Date",param)))$Date)
        # What date range do the measurements span?
          range(na.omit(subset(conc_df, conc_df$River==riv, select=c("Date",param)))$Date)
        # Plot # of mmnts per year
          obs_df <- na.omit(subset(conc_df, conc_df$River==riv, select=c("River","Date",param))) %>% group_by(year(Date)) %>% summarise(Count = n())
          names(obs_df) <- c("year","count")
          ggplot(obs_df, aes(y=count, x=year)) +
            # Plot
              geom_col() +
            # Limits
              scale_y_continuous(lim=c(0,24), breaks=seq(0,24,2)) +
              scale_x_continuous(lim=c(1970,2020), breaks=seq(1970,2020,10)) +
            # Themes
              theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    plot.background=element_rect(colour="white", size=1),
                    panel.border=element_rect(colour="black", fill=NA, size=1),
                    text=element_text(size=14)) +
              theme(axis.text.y=element_text(angle=0, hjust=0.5, colour="black")) +
              theme(axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0))) +
              theme(axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10))) +
              theme(axis.text.x=element_text(angle=0, hjust=0.5, colour="black")) +
              theme(axis.title.x=element_text(margin=margin(t=-10, r=0, b=0, l=0))) +  
              theme(plot.background=element_rect(fill='white')) +
              theme(panel.background=element_rect(fill='white')) +
              theme(legend.position="none") +
            # Labels
              labs(y=paste0("Number of observations (",riv," ",param,")"), x=NULL)
        # For which years is there at least one measurement?
          unique(year(na.omit(subset(conc_df, conc_df$River==riv, select=c("Date",param)))$Date))
        # Calculations
          tapply(X=conc_df$Hg, INDEX=conc_df$River, FUN=mean, na.rm=T)
          conc_df %>% group_by(River) %>% summarize(MeanHg = mean(Hg, na.rm=T)) # dyplr method for tapply
        