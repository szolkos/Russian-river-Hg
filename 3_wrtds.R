#===========================================================================================================#
# 3_wrtds.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | December 27, 2021
# Background: Weighted regressions on time, discharge, and season (WRTDS) modeling of conc and flux
#===========================================================================================================#

  # BACKGROUND ####

    # WRTDS information ####

      # WRTDS 
        # Provides flow-normalized (FN) concentration and flux estimates, in support of trend analysis (requires Q in m3/s) (Hirsch et al. 2010)
        # Developed for datasets with ≥ 200 samples, data spanning ≥ 10 years (Hirsch et al. 2015)
        # Can be run for datasets if ≥ 50 samples
      # WRTDS-K
        # Provides 'Kalman-filtered' best estimates of conc and flux which account for serial autocorrleation, in support of quantifying magnitudes (Zhang & Hirsch 2019)
      # Other
        # WRTDS-K and WRTDS will give similar results if sample frequency ~ ≥ 2 months; ~ ≤ 1 month may be different (Zhang & Hirsch 2019)

    # Parameters & Units ####
      # for conc and flux observations, and WRTDS results...
        # Q = discharge | m3/s
        # conc = concentration | mg/L
        # fn_conc = flow-normalized concentration | mg/L
        # flux = input flux | 10^6 kg/yr
        # fn_flux = flow-normalized flux | 10^6 kg/yr
        # gen_flux = Kalman flux (best estimate) | 10^6 kg/yr

  # FUNCTIONS ####
  
    # Compile and optionally export merged Q and C data for WRTDS modeling
      wrtds_data <- function(riv, param, export_yes_no, return_data_yes_no, print_summary_yes_no){
        
        # Set output directory
          #outdir <- paste0("/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Analyses/WRTDS/EGRET/", riv, "/", Sys.Date(), "/Input/")
          if(export_yes_no == "yes"){outdir <- paste0("/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Analyses/WRTDS/EGRET/",riv,"/")}
          if(export_yes_no == "no"){outdir <- ""}
  
        # Subset
          if(riv=="Onega" & param=="TSS" | riv=="Onega" & param=="Hg" | riv=="Onega" & param=="PHg"){year1 <- 1980; year2 <- 1992; wrtds_q <- on_q} # Same years for TSS, Hg, PHg; year2 determined by Q
          if(riv=="NorthDvina" & param=="TSS"){year1 <- 1979; year2 <- 2001; wrtds_q <- nd_q}
          if(riv=="NorthDvina" & param=="Hg" | riv=="NorthDvina" & param=="PHg"){year1 <- 1980; year2 <- 2001; wrtds_q <- nd_q} # Same years for Hg and PHg; year2 determined by Q
          if(riv=="Mezen" & param=="TSS"){year1 <- 1980; year2 <- 2018; wrtds_q <- me_q} #year2 determined by Q
          if(riv=="Mezen" & param=="Hg" | riv=="Mezen" & param=="PHg"){year1 <- 1980; year2 <- 2002; wrtds_q <- me_q} # year2 determined by Hg and PHg data availability
          if(riv=="Pechora" & param=="TSS"){year1 <- 1983; year2 <- 2018; wrtds_q <- pe_q}
          if(riv=="Pechora" & param=="Hg" | riv=="Pechora" & param=="PHg"){year1 <- 1983; year2 <- 2002; wrtds_q <- pe_q} 
          if(riv=="Ob" & param=="TSS"){year1 <- 1974; year2 <- 2008; wrtds_q <- ob_q}
          if(riv=="Ob" & param=="Hg" | riv=="Ob" & param=="PHg"){year1 <- 1979; year2 <- 2008; wrtds_q <- ob_q} # Same years for Hg and PHg
          if(riv=="Yenisey" & param=="TSS" | riv=="Yenisey" & param=="Hg" | riv=="Yenisey" & param=="PHg"){year1 <- 1980; year2 <- 2011; wrtds_q <- ye_q} # Same years for TSS, Hg, PHg; year1 determined by Q
          if(riv=="Lena" & param=="TSS"){year1 <- 2004; year2 <- 2018; wrtds_q <- le_q}
          if(riv=="Lena" & param=="Hg" | riv=="Lena" & param=="PHg"){year1 <- 2004; year2 <- 2011; wrtds_q <- le_q} # Same years for Hg and PHg
          if(riv=="Kolyma" & param=="TSS" | riv=="Kolyma" & param=="Hg" | riv=="Kolyma" & param=="PHg"){year1 <- 1982; year2 <- 1987; wrtds_q <- ko_q} # Same years for TSS, Hg, PHg; year1 determined by Q
          wrtds_riv <- merge(wrtds_q, subset(conc_df, conc_df$River==riv, select=c("Date","River","TSS","Hg","pHg")), all.x=T)
        # Q (daily, m3/s)
          wrtds_riv_q <- subset(wrtds_riv, select=c("Date","Q"), wrtds_riv$Year >= year1 & wrtds_riv$Year <= year2)
          names(wrtds_riv_q) <- c("date","Qdaily")
          wrtds_riv_q$date <- as.Date(wrtds_riv_q$date, format="%m/%d/%y")
          if(export_yes_no == "yes"){write.csv(wrtds_riv_q, paste0(dir, "Q_", riv, ".csv"), row.names=F)}
          if(print_summary_yes_no=="yes"){print(paste0("n daily discharge (m3/s) (", year1,"-", year2, "): ", length(na.omit(wrtds_riv_q$Qdaily)))); print(head(wrtds_riv_q))} # Display n
        # TSS (mg/L)
          if(param=="TSS"){
            wrtds_riv_tss <- na.omit(subset(wrtds_riv, select=c("Date","TSS"), wrtds_riv$Year >= year1 & wrtds_riv$Year <= year2))
            wrtds_riv_tss$rTSS <- ""
            names(wrtds_riv_tss) <- c("date","TSS","rTSS")
            wrtds_riv_tss <- na.omit(select(wrtds_riv_tss, date, rTSS, TSS))
            wrtds_riv_tss$date <- as.Date(wrtds_riv_tss$date, format="%m/%d/%y")
            if(export_yes_no == "yes"){write.csv(wrtds_riv_tss, paste0(dir, "S_", riv, "_TSS.csv"), row.names=F)}; 
            if(print_summary_yes_no=="yes"){print(paste0("n TSS (mg/L) (", year1,"-", year2, "): ", length(na.omit(wrtds_riv_tss$TSS)))); print(head(wrtds_riv_tss))}
            wrtds_riv_param <- wrtds_riv_tss
          }
        # Hg (µg/g)
          if(param=="Hg"){
            wrtds_riv_hg <- na.omit(subset(wrtds_riv, select=c("Date","Hg"), wrtds_riv$Year >= year1 & wrtds_riv$Year <= year2))
            wrtds_riv_hg$rHg <- ""
            names(wrtds_riv_hg) <- c("date","Hg","rHg")
            wrtds_riv_hg <- na.omit(select(wrtds_riv_hg, date, rHg, Hg))
            wrtds_riv_hg$date <- as.Date(wrtds_riv_hg$date, format="%m/%d/%y")
            if(export_yes_no == "yes"){write.csv(wrtds_riv_hg, paste0(dir, "S_", riv, "_Hg.csv"), row.names=F)}
            if(print_summary_yes_no=="yes"){print(paste0("n Hg (µg/g) (", year1,"-", year2, "): ", length(na.omit(wrtds_riv_hg$Hg)))); print(head(wrtds_riv_hg))}
            wrtds_riv_param <- wrtds_riv_hg
          }
        # PHg (ng/L)
          if(param=="PHg"){
            wrtds_riv_phg <- na.omit(subset(wrtds_riv, select=c("Date","pHg"), wrtds_riv$Year >= year1 & wrtds_riv$Year <= year2))
            wrtds_riv_phg$pHg <- as.numeric(format(wrtds_riv_phg$pHg/1000000, scientific=F, digits=10))
            wrtds_riv_phg$rPHg <- ""
            names(wrtds_riv_phg) <- c("date","PHg","rPHg")
            wrtds_riv_phg <- select(wrtds_riv_phg, date, rPHg, PHg)
            wrtds_riv_phg$date <- as.Date(wrtds_riv_phg$date, format="%m/%d/%y")
            if(export_yes_no == "yes"){write.csv(wrtds_riv_phg, paste0(dir, "S_", riv, "_PHg.csv"), row.names=F)}
            if(print_summary_yes_no=="yes"){print(paste0("n PHg (ng/L) (", year1,"-", year2, "): ", length(na.omit(wrtds_riv_phg$PHg)))); print(head(wrtds_riv_phg))}
            wrtds_riv_param <- wrtds_riv_phg
          }
          if(return_data_yes_no == "yes"){return(list(wrtds_riv_q, wrtds_riv_param))}
        }
    
    # Compile Hg and TSS conc and flux observations and WRTDS results
      wrtds_results <- function(riv){

        # Working directory
          wdir <- paste0("/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Analyses/WRTDS/EGRET/",riv,"/")
        # Compile data
          # TSS
            # Conc and fluxes
              tss_df <- read.csv(paste0(wdir, "annual_",riv,"_TSS.csv"), header=T, stringsAsFactors=F)
              names(tss_df) <- c("year","q","conc","fn_conc","flux","fn_flux","gen_flux")
            # Bootstrapped confidence intervals
              tss_ci_df <- read.csv(paste0(wdir, "CI_",riv,"_TSS.csv"), header=T, stringsAsFactors=F)
              names(tss_ci_df) <- c("year","fn_conc_lo","fn_conc_hi","fn_flux_lo","fn_flux_hi")
              tss_ci_df$year <- as.numeric(sub("\\..*", "", tss_ci_df$year))
            # Merge data
              tss_df <- merge(tss_df, tss_ci_df, by="year", all.x=T)
              tss_df <- subset(tss_df, select=c("year", "q", "conc", "fn_conc", "fn_conc_lo", "fn_conc_hi", "flux", "gen_flux", "fn_flux", "fn_flux_lo", "fn_flux_hi"))
              tss_df$param <- "TSS"
            # Unit conversions
              tss_df$flux <- tss_df$flux/1000 # Convert 10^6 kg/yr to Tg/yr
              tss_df$gen_flux <- tss_df$gen_flux/1000 # Convert 10^6 kg/yr to Tg/yr
              tss_df$fn_flux <- tss_df$fn_flux/1000 # Convert 10^6 kg/yr to Tg/yr
              tss_df$fn_flux_lo <- tss_df$fn_flux_lo*365/1000000000 # Convert kg/d to Tg/yr
              tss_df$fn_flux_hi <- tss_df$fn_flux_hi*365/1000000000 # Convert kg/d to Tg/yr
          # Hg  (turn off for Lena)
            # Conc and fluxes
              hg_df <- read.csv(paste0(wdir, "annual_",riv,"_Hg.csv"), header=T, stringsAsFactors=F)
              names(hg_df) <- c("year","q","conc","fn_conc","flux","fn_flux","gen_flux")
              hg_df <- droplevels(subset(hg_df, select=c("year","q","conc","fn_conc")))
            # Bootstrapped confidence intervals
              hg_ci_df <- read.csv(paste0(wdir, "CI_",riv,"_Hg.csv"), header=T, stringsAsFactors=F)
              names(hg_ci_df) <- c("year","fn_conc_lo","fn_conc_hi")
              hg_ci_df$year <- as.numeric(sub("\\..*", "", hg_ci_df$year))
            # Merge data
              hg_df <- merge(hg_df, hg_ci_df, by="year", all.x=T)
              hg_df <- subset(hg_df, select=c("year", "q", "conc", "fn_conc", "fn_conc_lo", "fn_conc_hi"))
              hg_df$param <- "Hg"
          # PHg fluxes  (turn off for Lena)
            # Conc and fluxes
              phg_df <- read.csv(paste0(wdir, "annual_",riv,"_PHg.csv"), header=T, stringsAsFactors=F)
              names(phg_df) <- c("year","q","conc","fn_conc","flux","fn_flux","gen_flux")
            # Bootstrapped confidence intervals
              hg_ci_df <- read.csv(paste0(wdir, "CI_",riv,"_PHg.csv"), header=T, stringsAsFactors=F)
              names(hg_ci_df) <- c("year","fn_conc_lo","fn_conc_hi","fn_flux_lo","fn_flux_hi")
              hg_ci_df$year <- as.numeric(sub("\\..*", "", hg_ci_df$year))
            # Merge data
              phg_df <- merge(phg_df, hg_ci_df, by="year", all.x=T)
              phg_df <- subset(phg_df, select=c("year", "q", "flux", "gen_flux", "fn_flux", "fn_flux_lo", "fn_flux_hi"))
              phg_df$param <- "PHg"
            # Unit conversions
              phg_df$flux <- phg_df$flux*1000 # Convert 10^6 kg/yr to Mg/yr (flux * 10^6 ÷ 1000, i.e., ÷ 1000)
              phg_df$gen_flux <- phg_df$gen_flux*1000 # Convert from 10^6 kg/yr to Mg/yr
              phg_df$fn_flux <- phg_df$fn_flux*1000 # Convert 10^6 kg/yr to Mg/yr
              phg_df$fn_flux_lo <- phg_df$fn_flux_lo*365/1000 # Convert kg/d to Mg/yr
              phg_df$fn_flux_hi <- phg_df$fn_flux_hi*365/1000 # Convert kg/d to Mg/yr
              
        # Return data
          return(list(tss_df, hg_df, phg_df))
      }
    
    # Plot Hg and TSS conc and flux observations and WRTDS results
      plot_wrtds_conc <- function(riv){
        
        # Data
          tss_df <- wrtds_results(riv=riv)[[1]] # TSS
          hg_df <- wrtds_results(riv=riv)[[2]] # Hg
        
        # Conditionals    
          if(riv=="Onega"){ylims_tss=c(0,40); ybreaks_tss=seq(0,40,10); ylims_phg=c(0,0.61); ybreaks_phg=seq(0,0.6,0.1); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          if(riv=="NorthDvina"){ylims_tss=c(0,50); ybreaks_tss=seq(0,50,10); ylims_phg=c(0,1.3); ybreaks_phg=seq(0,1.2,0.4); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          if(riv=="Mezen"){ylims_tss=c(0,50); ybreaks_tss=seq(0,50,10); ylims_phg=c(0,1); ybreaks_phg=seq(0,1,0.2); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          if(riv=="Pechora"){ylims_tss=c(0,60); ybreaks_tss=seq(0,60,20); ylims_phg=c(0,1.4); ybreaks_phg=seq(0,1.2,0.4); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          if(riv=="Ob"){ylims_tss=c(0,80); ybreaks_tss=seq(0,80,20); ylims_phg=c(0,0.82); ybreaks_phg=seq(0,0.8,0.2); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          if(riv=="Yenisey"){ylims_tss=c(0,30); ybreaks_tss=seq(0,30,5); ylims_phg=c(0,3); ybreaks_phg=seq(0,3,1); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          if(riv=="Lena"){ylims_tss=c(0,160); ybreaks_tss=seq(0,160,40); ylims_phg=c(0,1.4); ybreaks_phg=seq(0,1.4,0.4); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0; alpha_hg_ci_line=0}
          if(riv=="Kolyma"){ylims_tss=c(0,100); ybreaks_tss=seq(0,100,20); ylims_phg=c(0,1.3); ybreaks_phg=seq(0,1.2,0.4); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
          
        # Get mean measured concentrations and calculate scale factor for 2nd y-axis for TSS
          ann_conc_df <- ann_conc(riv)
          scale_factor <- max(na.omit(ann_conc_df$TSS_mean))/max(na.omit(ann_conc_df$Hg_mean))
          
        # Particulate mercury (Hg, ng/g) and total suspended solids (TSS, mg/L)
          plot_conc <- ggplot() +
            # Limits
              scale_y_continuous(lim=ylims_phg, breaks=ybreaks_phg, sec.axis=sec_axis(~.*scale_factor, breaks=ybreaks_tss, name=expression(TSS~(mg~L^-1)))) +
              scale_x_continuous(lim=c(1970,2020), breaks=seq(1970,2020,10)) +
            # Themes
              theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    plot.background=element_rect(colour="white", size=1),
                    panel.border=element_rect(colour="black", fill=NA, size=1),
                    text=element_text(size=22)) +
              theme(axis.text.y=element_text(angle=0, hjust=0.5, colour="black"), axis.text.y.right=element_text(color="grey40")) +
              theme(axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0))) +
              theme(axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10), color="grey40")) +
              theme(axis.text.x=element_text(angle=0, hjust=0.5, colour="black")) +
              theme(axis.title.x=element_text(margin=margin(t=-10, r=0, b=0, l=0))) +  
              theme(plot.background=element_rect(fill='white')) +
              theme(panel.background=element_rect(fill='white')) +
              theme(legend.position="none") +
            # TSS bootstrap CI from WRTDS
              geom_ribbon(data=tss_df, aes(x=year, ymin=fn_conc_lo/scale_factor, ymax=fn_conc_hi/scale_factor), fill="grey90", alpha=alpha_tss_ci_rib) +
              geom_line(data=tss_df, aes(x=year, y=fn_conc/scale_factor), color="white", size=1.2, alpha=alpha_tss_ci_line) +
            # TSS mean annual measured
              geom_errorbar(data=ann_conc_df, aes(x=Year, ymin=(TSS_mean-TSS_SE)/scale_factor, ymax=(TSS_mean+TSS_SE)/scale_factor), color="grey40", width=0, position=position_dodge()) +
              geom_point(data=ann_conc_df, aes(x=Year, y=TSS_mean/scale_factor), size=3.5, fill="white", color="grey40", shape=21, stroke=1.2) +
            # Hg bootstrap CI from WRTDS (turn off for Lena, also need to turn off Hg in 'wrtds_results')
              geom_ribbon(data=hg_df, aes(ymin=fn_conc_lo, ymax=fn_conc_hi, x=year), fill="grey10", alpha=alpha_hg_ci_rib) +
              geom_line(data=hg_df, aes(y=fn_conc, x=year), size=1.2, color="grey40", alpha=alpha_hg_ci_line) +
            # Hg mean annual measured
              geom_errorbar(data=ann_conc_df, aes(x=Year, ymin=Hg_mean-Hg_SE, ymax=Hg_mean+Hg_SE), color="black", width=0, position=position_dodge()) +
              geom_point(data=ann_conc_df, aes(x=Year, y=Hg_mean), size=3.5, fill="grey40", color="black", shape=21, stroke=1.2) +
            # Labels
              labs(y=PHg~(µg~g^-1), x=NULL)
          
            # Plot it
              print(plot_conc)
        }
              
      plot_wrtds_flux <- function(riv){
        
        # Run Lena and Kolyma manually for TSS, because no PHg trends analysis
          #riv <- "Lena"; scale_factor <- 3.57
          #riv <- "Kolyma"; scale_factor <- 5
        
        # Working directory
          wdir <- paste0("/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Analyses/WRTDS/EGRET/",riv,"/")
      
        # Compile data
          tss_df <- wrtds_results(riv=riv)[[1]]
          phg_df <- wrtds_results(riv=riv)[[3]]
          
        # Plotting
            
          # Conditionals
            if(riv=="Onega"){ylims_phg=c(0,0.25); ybreaks_phg=seq(0,0.25,0.05); ylims_tss=c(0,1); ybreaks_tss=seq(0,1,0.2); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
            if(riv=="NorthDvina"){ylims_phg=c(0,6); ybreaks_phg=seq(0,6,1); ylims_tss=c(0,5); ybreaks_tss=seq(0,5,1); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
            if(riv=="Mezen"){ylims_phg=c(0,0.62); ybreaks_phg=seq(0,0.6,0.1); ylims_tss=c(0,1); ybreaks_tss=seq(0,1,0.2); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
            if(riv=="Pechora"){ylims_tss=c(0,10); ybreaks_tss=seq(0,10,2); ylims_phg=c(0,4.3); ybreaks_phg=seq(0,4,1); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
            if(riv=="Ob"){ylims_phg=c(0,18); ybreaks_phg=seq(0,15,5); ylims_tss=c(0,50); ybreaks_tss=seq(0,50,10); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
            if(riv=="Yenisey"){ylims_phg=c(0,40); ybreaks_phg=seq(0,40,10); ylims_tss=c(0,12); ybreaks_tss=seq(0,12,4); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0.45; alpha_hg_ci_line=1}
            if(riv=="Lena"){ylims_phg=c(0,11); ybreaks_phg=seq(0,10,2); ylims_tss=c(0,40); ybreaks_tss=seq(0,40,10); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0; alpha_hg_ci_line=0}
            if(riv=="Kolyma"){ylims_phg=c(0,5); ybreaks_phg=seq(0,5,1); ylims_tss=c(0,25); ybreaks_tss=seq(0,25,5); alpha_tss_ci_rib=1; alpha_tss_ci_line=1; alpha_hg_ci_rib=0; alpha_hg_ci_line=0}
            
          # Get mean measured concentrations and calculate scale factor for 2nd y-axis for TSS
            
            scale_factor <- max(na.omit(tss_df$fn_flux))/max(na.omit(phg_df$fn_flux))
      
          # Conditionals and data points from literature
            # Null data
              if(riv=="Onega"){x1=0; y1=0; x2=0; y2=0; x3=0; y3=0; x4=0; y4=0; se1=0; se2=0; se3=0; se4=0; rng=0; zea_cent=0; zea_rng=0; sea_cent=0; sea_rng=0}
              if(riv=="NorthDvina"){x1=0; y1=0; x2=0; y2=0; x3=0; y3=0; x4=0; y4=0; se1=0; se2=0; se3=0; se4=0; rng=0; zea_cent=0; zea_rng=0; sea_cent=0; sea_rng=0}
              if(riv=="Mezen"){x1=0; y1=0; x2=0; y2=0; x3=0; y3=0; x4=0; y4=0; se1=0; se2=0; se3=0; rng=0; se4=0; zea_cent=0; zea_rng=0; sea_cent=0; sea_rng=0}
              if(riv=="Pechora"){x1=0; y1=0; x2=0; y2=0; x3=0; y3=0; x4=0; y4=0; se1=0; se2=0; se3=0; rng=0; se4=0; zea_cent=0; zea_rng=0; sea_cent=0; sea_rng=0}
            # Zolkos et al. 2018 ES&T center timepoint
              zea_cent <- 2012+(2017-2012)/2; zea_rng <- zea_cent-2012
            # Ob'
              if(riv=="Ob"){
                # THg- Zolkos et al. 2020; SE = 647/√34
                  zea_cent <- zea_cent; zea_rng <- zea_rng; x1=zea_cent; y1=2421/1000; rng=zea_rng; se1=((647/1000)/sqrt(34))/2
                # TSS- Zolkos et al. 2020 (mean ± SD); NOTE: plots to primary axis, so need to scale TSS values w/ "scale_factor", from above
                  zea_cent_tss <- 2012+(2017-2012)/2; zea_rng_tss <- zea_cent_tss-2012; x2=zea_cent_tss; y2=7.9/scale_factor; rng=zea_rng_tss; se2=2.3/scale_factor
                # THg- Sonke et al. 2018 PNAS
                  sea_cent <- 2012+(2016-2012)/2; sea_rng <- sea_cent-2012; x3=sea_cent ; y3=1586/1000; se3=0
                # THg- Lim et al. 2019 Env. Poll.
                  x4=2016; y4=1.5; se4=0.15
              } 
            # Yenisey
              if(riv=="Yenisey"){
                # Zolkos et al. 2020; SE = 567/√33
                  x1=zea_cent; y1=3642/1000; rng=zea_rng; se1=((567/1000)/sqrt(33))/2
                # TSS- Zolkos et al. 2020 (mean ± SD)
                  zea_cent_tss <- 2012+(2017-2012)/2; zea_rng_tss <- zea_cent_tss-2012; x2=zea_cent_tss; y2=2.3/scale_factor; rng=zea_rng_tss; se2=0.6/scale_factor
                # Sonke et al. 2018 PNAS
                  sea_cent <- 2012+(2016-2012)/2; sea_rng <- sea_cent-2012; x3=sea_cent ; y3=2867/1000 ; se3=0
                # Null data
                  x4=0; y4=0; se4=0
              }
            # Lena
              if(riv=="Lena"){
                # Zolkos et al. 2020; SE = 406/√33
                  x1=zea_cent; y1=6591/1000; rng=zea_rng; se1=((406/1000)/sqrt(33))/2
                # TSS- Zolkos et al. 2020 (mean ± SD)
                  zea_cent_tss <- 2012+(2017-2012)/2; zea_rng_tss <- zea_cent_tss-2012; x2=zea_cent_tss; y2=23.1/scale_factor; rng=zea_rng_tss; se2=2.0/scale_factor
                # Sonke et al. 2018 PNAS
                  sea_cent <- 2012+(2016-2012)/2; sea_rng <- sea_cent-2012; x3=sea_cent ; y3=2430/1000 ; se3=0
                # Null data
                  x4=0; y4=0; se4=0
              }
            # Kolyma
              if(riv=="Kolyma"){
                # Zolkos et al. 2020; SE = 565/√29
                  x1=zea_cent; y1=1107/1000; rng=zea_rng; se1=((565/1000)/sqrt(29))/2
                # TSS- Zolkos et al. 2020 (mean ± SD)
                  zea_cent_tss <- 2012+(2017-2012)/2; zea_rng_tss <- zea_cent_tss-2012; x2=zea_cent_tss; y2=12.2/scale_factor; rng=zea_rng_tss; se2=11.1/scale_factor
                # Sonke et al. 2018 PNAS
                  sea_cent <- 2012+(2016-2012)/2; sea_rng <- sea_cent-2012; x3=sea_cent ; y3=483/1000 ; se3=0
                # Null data
                  x4=0; y4=0; se4=0
              }
            
        # Particulate mercury (Hg, ng/g) and total suspended solids (TSS, mg/L)
          plot_flux <- ggplot() +
            
            # Limits
              #scale_y_continuous(lim=ylims_phg, breaks=ybreaks_phg) +
              scale_y_continuous(lim=ylims_phg, breaks=ybreaks_phg, sec.axis=sec_axis(~.*scale_factor, breaks=ybreaks_tss, name=expression(TSS~(Tg~y^-1)))) +
              scale_x_continuous(lim=c(1970,2020), breaks=seq(1970,2020,10)) +
            
            # Themes
              theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    plot.background=element_rect(colour="white", size=1),
                    panel.border=element_rect(colour="black", fill=NA, size=1),
                    text=element_text(size=22)) +
              theme(axis.text.y=element_text(angle=0, hjust=0.5, colour="black"), axis.text.y.right=element_text(color="grey40")) +
              theme(axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0))) +
              theme(axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10), color="grey40")) +
              theme(axis.text.x=element_text(angle=0, hjust=0.5, colour="black")) +
              theme(axis.title.x=element_text(margin=margin(t=-10, r=0, b=0, l=0))) +  
              theme(plot.background=element_rect(fill='white')) +
              theme(panel.background=element_rect(fill='white')) +
              theme(legend.position="none") +
              
            # WRTDS
              # TSS bootstrap CI from WRTDS
                geom_ribbon(data=tss_df, aes(x=year, ymin=fn_flux_lo/scale_factor, ymax=fn_flux_hi/scale_factor), fill="grey90", alpha=alpha_tss_ci_rib) +
                geom_line(data=tss_df, aes(x=year, y=fn_flux/scale_factor), color="white", size=1.2, alpha=alpha_tss_ci_line) +
              # TSS mean annual estimated from WRTDS (Gen [Kalman] fluxes)
                geom_point(data=tss_df, aes(x=year, y=gen_flux/scale_factor), size=3.5, fill="white", color="grey40", shape=21, stroke=1.2) +
              # Hg bootstrap CI from WRTDS - (turn off for Lena & Kolyma)
                geom_ribbon(data=phg_df, aes(ymin=fn_flux_lo, ymax=fn_flux_hi, x=year), fill="grey10", alpha=alpha_hg_ci_rib) +
                geom_line(data=phg_df, aes(y=fn_flux, x=year), size=1.2, color="grey40", alpha=1) +
              # Hg mean annual estimated from WRTDS (Gen [Kalman] fluxes) - (turn off for Lena & Kolyma)
                geom_point(data=phg_df, aes(x=year, y=gen_flux), size=3.5, fill="grey40", color="black", shape=21, stroke=1.2) +
        
            # Labels
              labs(y=PHg~(Mg~y^-1), x=NULL) +
            
            # Data points from literature
              if(riv=="Ob" | riv=="Yenisey" | riv=="Lena" | riv=="Kolyma"){
                # THg- Zolkos et al. 2018 ES&T
                  list(geom_errorbar(data=tss_df, aes(x=x1, y=y1, ymin=y1-se1, ymax=y1+se1), color="black", width=0, position=position_dodge()),
                  geom_errorbarh(data=tss_df, aes(x=x1, y=y1, xmin=x1-zea_rng, xmax=x1+zea_rng), color="black", height=0, position=position_dodge()),
                  geom_point(data=tss_df, aes(x=x1, y=y1, size=1.2), fill="grey20", color="black", shape=23, size=2.8, stroke=1),
                # TSS- Zolkos et al. 2018 ES&T
                  geom_errorbar(data=tss_df, aes(x=x2, y=y2, ymin=y2-se2, ymax=y2+se2), color="black", width=0, position=position_dodge()),
                  geom_errorbarh(data=tss_df, aes(x=x2, y=y2, xmin=x2-zea_rng_tss, xmax=x2+zea_rng_tss), color="black", height=0, position=position_dodge()),
                  geom_point(data=tss_df, aes(x=x2, y=y2, size=1.2), fill="white", color="black", shape=23, size=2.8, stroke=1),
                # THg- Sonke et al. 2018 PNAS
                  geom_errorbarh(data=tss_df, aes(x=x3, y=y3, xmin=x3-sea_rng, xmax=x3+sea_rng), color="black", height=0, position=position_dodge()),
                  geom_point(data=tss_df, aes(x=x3, y=y3, size=1.2), fill="grey20", color="black", shape=24, size=2.8, stroke=1),
                # THg- Lim et al. 2019 Env. Poll.
                  geom_errorbar(data=tss_df, aes(x=x4, ymin=y4-se4, ymax=y4+se4), color="black", width=0, position=position_dodge()),
                  geom_point(data=tss_df, aes(x=x4, y=y4, size=1.2), fill="grey20", color="black", shape=22, size=2.8, stroke=1))}
          
          # Plot it
            print(plot_flux)
            
      }
    
    # Peak discharge magnitude & trends
      Q_compile <- function(riv, param, plot_yn, plot_save_yn){

        # Subset Q data
          qdf <- wrtds_data(riv=riv, param=param, export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="no")[[1]]
          head(qdf); tail(qdf); dim(qdf) # Inspect
          qdf$year <- year(qdf$date) # Add year to df
          
        # Table with total annual discharge (km3) and annual discharge peak magnitude (km3)
          years <- unique(qdf$year)
          df <- as.data.frame(matrix(ncol=4, nrow=length(years)))
          names(df) <- c("river","year","q_tot_ann","q_peak_mag")
          j=1
          for(i in years){
            # Subset data
              qdf_sub <- droplevels(subset(qdf, qdf$year==i)) # Subset data for year i
              pcntl_95 <- quantile(qdf_sub$Qdaily, 0.95) # Value of 95th percentile
              qdf_pcntl_95 <- droplevels(subset(qdf_sub, qdf_sub$Qdaily >= pcntl_95)) # Subset Q values within 95th percentile
            # Store to df
              df[j,1] <- riv
              df[j,2] <- i
              df[j,3] <- sum((as.numeric(qdf_sub$Qdaily)*86400)/1e9) # Total annual discharge, converted from m3 to km3
              df[j,4] <- sum((as.numeric(qdf_pcntl_95$Qdaily)*86400)/1e9) # Sum Q values within 95th percentile
              j <- j+1
          }
          
        # Plot

          # Calculate limits
            q_tot_ann_min <- round(min(df$q_tot_ann)-(min(df$q_tot_ann)*0.1),0)
            q_tot_ann_max <- round(max(df$q_tot_ann)+(max(df$q_tot_ann)*0.1),0)
          
          plot_q_tot_ann <- ggplot() +
            # Limits
              scale_y_continuous(lim=c(q_tot_ann_min,q_tot_ann_max), labels=NotFancy) +
              scale_x_continuous(lim=c(1970,2020), breaks=seq(1970,2020,10)) +
            # Themes
              theme(plot.margin=unit(c(0.15,0.12,0.12,0.12), "in"),
                    plot.background=element_rect(colour="white", fill="white", size=1),
                    panel.background=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.border=element_rect(colour="black", fill=NA, size=1)) +
              theme(axis.text.y=element_text(size=13, angle=0, hjust=0.5, colour="black"),
                    axis.text.x=element_text(size=13, angle=0, hjust=0.5, colour="black"),
                    axis.title.y.left=element_text(size=15, margin=margin(t=0, r=10, b=0, l=0)),
                    axis.title.x=element_blank(),
                    plot.title=element_text(size=13, angle=0, hjust=0, face="bold", colour="black")) +
              theme(legend.position="none") +
            # Data
              geom_point(data=df, aes(x=year, y=q_tot_ann), size=3.5, fill="grey40", color="black", shape=21, stroke=1.2) +
            # Labels
              labs(title=riv, y=Total~annual~italic(Q)~(km^3))
          
          # Calculate limits
            q_peak_mag_min <- round(min(df$q_peak_mag)-(min(df$q_peak_mag)*0.1),0)
            q_peak_mag_max <- round(max(df$q_peak_mag)+(max(df$q_peak_mag)*0.1),0)
            
          plot_q_peak_mag <- ggplot() +
            # Limits
              scale_y_continuous(lim=c(q_peak_mag_min,q_peak_mag_max), labels=NotFancy) +
              scale_x_continuous(lim=c(1970,2020), breaks=seq(1970,2020,10)) +
            # Themes
              theme(plot.margin=unit(c(0.15,0.12,0.12,0.12), "in"),
                    plot.background=element_rect(colour="white", fill="white", size=1),
                    panel.background=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.border=element_rect(colour="black", fill=NA, size=1)) +
              theme(axis.text.y=element_text(size=13, angle=0, hjust=0.5, colour="black"),
                    axis.text.x=element_text(size=13, angle=0, hjust=0.5, colour="black"),
                    axis.title.y.left=element_text(size=15, margin=margin(t=0, r=10, b=0, l=0)),
                    axis.title.x=element_blank(),
                    plot.title=element_text(size=13, angle=0, hjust=0, face="bold", colour="black")) +
              theme(legend.position="none") +
            # Data
              geom_point(data=df, aes(x=year, y=q_peak_mag), size=3.5, fill="grey40", color="black", shape=21, stroke=1.2) +
            # Labels
              labs(y=Peak~italic(Q)~magnitude~(km^3))
              
          # Plot it
            q_grob <- arrangeGrob(plot_q_tot_ann, plot_q_peak_mag, ncol=1)
            if(plot_yn=="yes"){plot(q_grob)}
            if(plot_save_yn=="yes"){
              out <- "/Users/szolkos/Documents/Research/Projects/Russian Arctic river TSS Hg/Manuscript drafts/2022-01-30-rvsd/Figures/Q trends/"
              ggsave(filename=paste0(out,riv,"_",Sys.Date(),".pdf"), plot=q_grob, height=5.5, width=6)
            }
          
        # Return data
          return(list(df))
      }
      
      Q_trends <- function(param, export_yn){

        # Table with trends in total annual discharge and discharge peak magnitude
          df <- as.data.frame(matrix(ncol=7, nrow=length(rus_rivers)))
          names(df) <- c("river","q_tot_ann_change_int","q_tot_ann_change_mag","q_tot_ann_sig","q_peak_mag_int","q_peak_mag_change_mag","q_peak_mag_sig")
          j=1
          for(i in rus_rivers){
            # Q data
              df_sub <- as.data.frame(Q_compile(i, param, plot_yn="no", plot_save_yn="no"))
            # Trends
              df[j,1] <- i
              # zyp.trend.vector outputs: 2 = trend (annual change), 11 = intercept, 3 = trendp (total change), 5 = tau, 6 = p
              df[j,2] <- zyp.trend.vector(y=as.numeric(as.matrix(subset(df_sub, select="q_tot_ann"))), method="yuepilon", conf.intervals=T)[11]
              df[j,3] <- zyp.trend.vector(y=as.numeric(as.matrix(subset(df_sub, select="q_tot_ann"))), method="yuepilon", conf.intervals=T)[3]
              df[j,4] <- zyp.trend.vector(y=as.numeric(as.matrix(subset(df_sub, select="q_tot_ann"))), method="yuepilon", conf.intervals=T)[6]
              df[j,5] <- zyp.trend.vector(y=as.numeric(as.matrix(subset(df_sub, select="q_peak_mag"))), method="yuepilon", conf.intervals=T)[11]
              df[j,6] <- zyp.trend.vector(y=as.numeric(as.matrix(subset(df_sub, select="q_peak_mag"))), method="yuepilon", conf.intervals=T)[3]
              df[j,7] <- zyp.trend.vector(y=as.numeric(as.matrix(subset(df_sub, select="q_peak_mag"))), method="yuepilon", conf.intervals=T)[6]
            j <- j+1
          }
          
        # Export
          if(export_yn=="yes"){
            write.csv(df, paste0(dir,"Q_trends_",param,"_",Sys.Date(),".csv"), row.names=F)
          }
          
        # Return data
          return(df)
      }

      Q_t1_t2 <- function(param){
        Names <- c("river","q_tot_ann_t1","q_tot_ann_t2","q_peak_mag_t1","q_peak_mag_t2")
        df <- as.data.frame(matrix(ncol=length(Names), nrow=length(rus_rivers)))
        names(df) <- Names
        j=1
        for(i in rus_rivers){
          df[j,1] <- i
          df[j,2] <- first(Q_compile(i, param, plot_yn="no", plot_save_yn="no")[[1]]$q_tot_ann)
          df[j,3] <- last(Q_compile(i, param, plot_yn="no", plot_save_yn="no")[[1]]$q_tot_ann)
          df[j,4] <- first(Q_compile(i, param, plot_yn="no", plot_save_yn="no")[[1]]$q_peak_mag)
          df[j,5] <- last(Q_compile(i, param, plot_yn="no", plot_save_yn="no")[[1]]$q_peak_mag)
          j <- j+1
        }
        return(df)
      }
    
    # Yields
      yields <- function(riv, param){
        
        # Calculate TSS and PHg yields using GEN fluxes
          # Store watershed areas
            if(riv=="Onega"){shed_area <- 55700}
            if(riv=="NorthDvina"){shed_area <- 348000}
            if(riv=="Mezen"){shed_area <- 564000}
            if(riv=="Pechora"){shed_area <- 248000}
            if(riv=="Ob"){shed_area <- 2990000}
            if(riv=="Yenisey"){shed_area <- 2400000}
            if(riv=="Lena"){shed_area <- 2430000}
            if(riv=="Kolyma"){shed_area <- 530000}
          # Compile data
            if(param=="TSS"){
              ylds <- wrtds_results(riv=riv)[[1]]
              ylds$shed_area <- shed_area
              ylds$yield <- ylds$gen_flux*1000000/ylds$shed_area # *1000000 converts TSS fluxes from Tg to Mg
              ylds$river <- riv
              ylds <- droplevels(subset(ylds, select=c("river", "year", "yield")))
              names(ylds) <- c("river", "year", "tss_yld")
            }
            if(param=="PHg"){
              ylds <- wrtds_results(riv=riv)[[3]]
              ylds$shed_area <- shed_area
              ylds$yield <- ylds$gen_flux*1000000/ylds$shed_area # *1000000 converts PHg fluxes from Mg to g
              ylds$river <- riv
              ylds <- droplevels(subset(ylds, select=c("river", "year", "yield")))
              names(ylds) <- c("river", "year", "phg_yld")
            }
        # Return data
          return(ylds)
          
      }
    
    # Compile data for Fig. S1 (TSS and PHg sampling days as points superimposed on hydrographs of all 8 rivers)
      Fig_S1_data <- function(riv){
        # Compile dates of Q observation and parameter sampling
          dates_q <- wrtds_data(riv=riv, param="TSS", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="no")[[1]]
          dates_tss <- wrtds_data(riv=riv, param="TSS", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="no")[[2]]
          dates_phg <- wrtds_data(riv=riv, param="PHg", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="no")[[2]]
        # Combine into df
          dates_df <- merge(dates_q, dates_tss, by="date", all.x=T)
          dates_df <- merge(dates_df, dates_phg, by="date", all.x=T)
          dates_df$river <- riv
          dates_df <- droplevels(subset(dates_df, select=c("river","date","Qdaily","TSS","PHg")))
          names(dates_df) <- c("river","date","qdaily","tss","phg")
          dates_df$qdaily <- dates_df$qdaily*86400
        # Return
          return(dates_df)
      }
      
        
  # DATA PROCESSING AND ANALYSIS ####
        
    # Compile and Export data for WRTDS ####
  
      # Select river & parameter
        fig_S1_df <- rbind(Fig_S1_data("Onega"), Fig_S1_data("NorthDvina"), Fig_S1_data("Mezen"), Fig_S1_data("Pechora"), Fig_S1_data("Ob"), Fig_S1_data("Yenisey"), Fig_S1_data("Lena"), Fig_S1_data("Kolyma"))
        fig_S1_df$river <- factor(fig_S1_df$river, levels=rus_rivers)
      # Optionally export
        #write.csv(fig_S1_df, paste0(dir,"fig_S1_df.csv"), row.names=F)
        
    # Compile and plot results from WRTDS ####
        
      # Select river and date on which WRTDS outputs were generated, then plot & save conc & flux
        riv <- "Yenisey"
      # Conc
        plot_wrtds_conc(riv=riv)
        ggsave(paste0(dir,"_",Sys.Date(),".pdf"), height=4, width=8)
      # Flux
        plot_wrtds_flux(riv=riv)
        ggsave(paste0(dir,"_",Sys.Date(),".pdf"), height=4, width=8)
        
        
    # Statistics ####
        
      # Summary (n obs) for each parameter used in WRTDS modeling, by river
        wrtds_data(riv="Kolyma", param="PHg", export_yes_no="no", return_data_yes_no="no", print_summary_yes_no="no")
        
      # Start & end dates of sample collection for WRTDS modeling durations for TSS and PHg
        riv <- "Kolyma"
        first(wrtds_data(riv=riv, param="TSS", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="")[[2]])
        last(wrtds_data(riv=riv, param="TSS", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="")[[2]])
        first(wrtds_data(riv=riv, param="Hg", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="")[[2]])
        last(wrtds_data(riv=riv, param="Hg", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="")[[2]])
        first(wrtds_data(riv=riv, param="TSS", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="")[[1]])
        last(wrtds_data(riv=riv, param="TSS", export_yes_no="no", return_data_yes_no="yes", print_summary_yes_no="")[[1]])
      
      # Summary statistics on Q and concentration database
        riv <- "Kolyma"
        param <- "TSS"
        dim(na.omit(subset(conc_df, conc_df$River==riv, select=c("Date","River",param))))
        head(na.omit(subset(conc_df, conc_df$River==riv, select=c("Date","River",param))))[1,]
        tail(na.omit(subset(conc_df, conc_df$River==riv, select=c("Date","River",param))))[6,]
        
        dim(na.omit(subset(wrtds_df, wrtds_df$River==riv, select=c("Date","River",param))))
        head(na.omit(subset(wrtds_df, wrtds_df$River==riv, select=c("Date","River",param))))
        tail(na.omit(subset(wrtds_df, wrtds_df$River==riv, select=c("Date","River",param))))
        
      # Determine peak discharge magnitude & trends, optionally plot (run for param=TSS, which covers the full discharge records for both TSS and Hg)
        Q_compile("Onega", "TSS", plot_yn="no", plot_save_yn="no")[[1]]
        Q_trends("PHg", export_yn="no")
          
      # Total annual Q and peak Q magnitude for t1 and t2 for each river
        param <- "TSS"
        Q_t1_t2_df <- Q_t1_t2(param)
      # Optionally export
        #write.csv(Q_t1_t2_df, paste0(dir,"Q_t1_t2_df_",param,"_",Sys.Date(),".csv"), row.names=F)
        
      # Compare yields among rivers
        # Compile data
          yld_tss <- rbind(yields("Onega","TSS"), yields("NorthDvina","TSS"), yields("Mezen","TSS"), yields("Pechora","TSS"), yields("Ob","TSS"), yields("Yenisey","TSS"), yields("Lena","TSS"), yields("Kolyma","TSS"))
          yld_tss <- droplevels(subset(yld_tss, select=c("river","year","tss_yld")))
          yld_phg <- rbind(yields("Onega","PHg"), yields("NorthDvina","PHg"), yields("Mezen","PHg"), yields("Pechora","PHg"), yields("Ob","PHg"), yields("Yenisey","PHg"), yields("Lena","PHg"), yields("Kolyma","PHg"))
          yld_phg <- droplevels(subset(yld_phg, select=c("river","year","phg_yld"), yld_phg$river!="Lena" & yld_phg$river!="Kolyma"))
        # Merge data
          ylds <- merge(yld_tss, yld_phg, by=c("river","year"), all.x=T)
          ylds$river <- factor(ylds$river, levels=rus_rivers)
          
      # Compare means among groups (i.e. rivers) 
        # Calculate means and SE by river
          tapply(ylds$tss_yld, ylds$river, FUN=mean, na.rm=T); tapply(ylds$tss_yld, ylds$river, FUN=std.error, na.rm=T)
          tapply(ylds$phg_yld, ylds$river, FUN=mean, na.rm=T); tapply(ylds$phg_yld, ylds$river, FUN=std.error, na.rm=T)
        # Test data for assumptions of ANOVA: (i) equal variance, (ii) normally distributed, (iii) independent residuals
          bartlett.test( log(tss_yld)~river, ylds) # Variances unequal, even w/ ln transformation
          bartlett.test( log(phg_yld)~river, ylds) # Variances unequal, even w/ ln transformation
        # Welch's ANOVA for unequal variances | http://www.biostathandbook.com/onewayanova.html#welch, https://www.rdocumentation.org/packages/onewaytests/versions/2.4/topics/welch.test
          # TSS
            tss_for_welch <- droplevels(subset(ylds, select=c("river","year","tss_yld"))) # subset data
            tss_for_welch$ln_tss_yld <- log(tss_for_welch$tss_yld) # add ln-transformed TSS yield to df
            nor.test(ln_tss_yld ~ river, data=tss_for_welch, method="SW", plot="qqplot-histogram") # test for normality in residuals, ln transform data as needed
            ylds_mdl_tss <- welch.test(ln_tss_yld ~ river, tss_for_welch)
            paircomp(ylds_mdl_tss, adjust.method="bonferroni")
          # PHg
            phg_for_welch <- droplevels(na.omit(subset(ylds, select=c("river","year","phg_yld"))))
            phg_for_welch$ln_phg_yld <- log(phg_for_welch$phg_yld)
            nor.test(ln_phg_yld ~ river, data=phg_for_welch, method="SW", plot="qqplot-histogram")
            ylds_mdl_phg <- welch.test(ln_phg_yld ~ river, phg_for_welch)
            paircomp(ylds_mdl_phg, adjust.method="bonferroni") # pairwise comparisons
        # Mean TSS yields in western vs. eastern watersheds
          mean(na.omit(subset(ylds, select=c("year","river","tss_yld"), ylds$river != "Lena" | ylds$river != "Kolyma"))$tss_yld) # Western (mean)
          std.error (na.omit(subset(ylds, select=c("year","river","tss_yld"), ylds$river != "Lena" | ylds$river != "Kolyma"))$tss_yld) # Western (SE)
          mean(subset(ylds, select=c("year","river","tss_yld"), ylds$river == "Lena" | ylds$river == "Kolyma")$tss_yld) # Eastern (mean)
          std.error(subset(ylds, select=c("year","river","tss_yld"), ylds$river == "Lena" | ylds$river == "Kolyma")$tss_yld) # Eastern (SE)
            