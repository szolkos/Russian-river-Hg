#===========================================================================================================#
# 7_drm.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) | December 22, 2021
# Background: Dynamic regression model (w/ lags) for predicting [PHg] among Russian rivers
#===========================================================================================================#

  # BACKGROUND ####

    # Parameters & Units ####
  
      # river = Onega, NorthDvina, Mezen, Pechora, Ob, Yenisey, Lena, Kolyma
      # year = sampling year
      # phg = mean annual particulate mercury (Hg) concentration (µg/g)
      # phg_fn = particulate mercury (flow-normalized)
      # tss = mean annual total suspended solids concentration (mg/L)
      # tss_fn = total suspended solids (flow-normalized)
      # maat = mean annual air temperature (˚C)
      # tap = total annual precipitation (m)
      # regional = annual Hg deposition from regional (EU + USSR) primary anthropogenic Hg emissions
      # regional_n = regional, lagged by n years
      # global = annual Hg deposition from global primary anthropogenic Hg emissions
      # global_n = global, lagged by n years


  # FUNCTIONS ####
  
    # Exclude columns using 'subset()'
      `%ni%` <- Negate(`%in%`)
    
    # Compile and return flow-normalized concentration or flux for desired river
      fn_compile <- function(riv, conc_or_flux){
        
        # Hg conc and fluxes
          hg_df <- read.csv(paste0(dir, "WRTDS_EGRET/annual_",riv,"_Hg.csv"), header=T, stringsAsFactors=F)
          names(hg_df) <- c("year","q","conc","fn_conc_hg","flux","fn_flux_hg","gen_flux")
          hg_df <- select(hg_df, year, fn_conc_hg, fn_flux_hg)
        # TSS conc and fluxes
          tss_df <- read.csv(paste0(dir, "WRTDS_EGRET/annual_",riv,"_TSS.csv"), header=T, stringsAsFactors=F)
          names(tss_df) <- c("year","q","conc","fn_conc_tss","flux","fn_flux_tss","gen_flux")
          tss_df <- select(tss_df, year, fn_conc_tss, fn_flux_tss)
        # Merge Hg and TSS conc and flux data
          fn_df <- merge(hg_df, tss_df, by="year", all.x=T)
          fn_df$river <- riv
          fn_conc <- droplevels(select(fn_df, year, river, fn_conc_hg, fn_conc_tss))
          fn_flux <- droplevels(select(fn_df, year, river, fn_flux_hg, fn_flux_tss))
        # Return data
          if(conc_or_flux=="conc"){return(fn_conc)}
          if(conc_or_flux=="flux"){return(fn_flux)}
          
      }
    
      
  # DATA PROCESSING & ANALYSIS ####
    
    # Extract flow-normalized (FN) Hg and TSS conc data
      on_fn <- fn_compile("Onega", "conc")
      nd_fn <- fn_compile("NorthDvina", "conc")
      me_fn <- fn_compile("Mezen", "conc")
      pe_fn <- fn_compile("Pechora", "conc")
      ob_fn <- fn_compile("Ob", "conc")
      ye_fn <- fn_compile("Yenisey", "conc")
      #le_fn <- fn_compile("Lena", "conc")
      #ko_fn <- fn_compile("Kolyma", "conc")
    # Compile FN Hg and TSS conc data
      fn_df <- rbind(on_fn, nd_fn, me_fn, pe_fn, ob_fn, ye_fn)
      unique(fn_df$river)
      fn_df$river <- factor(fn_df$river, levels=c("Onega","NorthDvina","Mezen","Pechora","Ob","Yenisey"))
      summary(fn_df)
    # Merge flow-normalized and climate data
      era_sub <- droplevels(subset(era, select=c("year","river","maat","tap")))
      drm_df <- merge(fn_df, era_sub, by=c("year","river"), all.x=T)
    # Merge climate and reservoir data
      drm_df <- merge(drm_df, res_cum_cap, by=c("river","year"), all.x=T)
    # Merge Hg deposition data
      Hg_dep_sub <- droplevels(subset(Hg_dep, select=c("river","year","regional_0","regional_5","regional_10","regional_15","global_0","global_5","global_10","global_15")))
      drm_df <- merge(drm_df, Hg_dep_sub, by=c("river","year"), all.x=T)
    # Set factor levels
      drm_df$river <- factor(drm_df$river, levels=c("Onega", "NorthDvina", "Mezen", "Pechora", "Ob", "Yenisey"))
               
      
    # Dynamic Regression Model ####
      
      # Note: model development is hands-on & iterative, thus the following steps are done without the aid of homemade functions
        
      # Subset data
        riv <- "Yenisey"
        arima_pre <- droplevels(drm_df %>% filter(river == riv) %>% select(fn_conc_hg, year, fn_conc_tss, maat, tap, res_cum_cap_20, res_cum_cap_25, regional_0, regional_5, regional_10, regional_15, global_0, global_5, global_10, global_15))
        
      # Examine covariate relationships
        #ggpairs(arima_pre)
        
      # Durbin-Watson test of autocorrelation in dependent variable; p < 0.05 indicates autocorrelation
        durbinWatsonTest(lm(arima_pre$fn_conc_hg ~ arima_pre$year)) # Test residuals
        
      # KPSS informs n differences (e.g., 1st-order AR, i.e. AR(1)) | https://otexts.com/fpp3/non-seasonal-arima.html#fig:egyptacf
        # Estimation and order selection (i.e., determination of p,d,q) | https://otexts.com/fpp3/non-seasonal-arima.html, https://otexts.com/fpp3/arima-estimation.html
          # p = order of the autoregressive part; d = degree of first-order differencing; q = order of the moving average part
            # pdq(p=1:n, d=1:n, q=0:n) determines optimal p,d,q based on AICc, but user should first determine appropriate range of values for d,q
      # Use KPSS to test covariates for stationarity (https://otexts.com/fpp3/stationarity.html | https://otexts.com/fpp3/estimation.html)
        # Unit root test: Ho: data are stationary; if p < 0.05, then Ho rejected and data are non-stationary and thus differencing may be required
          unitroot_kpss(as.matrix(arima_pre %>% select(fn_conc_hg)))
        # If not stationary (p < 0.05), determine appropriate number of first differences
          ndiffs <- unitroot_ndiffs(as.matrix(arima_pre %>% select(fn_conc_hg))); print(ndiffs)
        # From previous, apply n 1st differences, re-test stationarity
          unitroot_kpss(diff(as.matrix(arima_pre %>% select(fn_conc_hg)), ndiffs))
        
      # Create ts series tibble (tsibble) of data, to assess y (i.e., [PHg]) using ACF & PACF
        arima_pre <- tsibble(year=arima_pre$year, fn_conc_hg=arima_pre$fn_conc_hg, fn_conc_tss=arima_pre$fn_conc_tss, maat=arima_pre$maat, tap=arima_pre$tap, res_cum_cap_20=arima_pre$res_cum_cap_20, res_cum_cap_25=arima_pre$res_cum_cap_25, regional_0=arima_pre$regional_0, regional_5=arima_pre$regional_5, regional_10=arima_pre$regional_10, regional_15=arima_pre$regional_15, global_0=arima_pre$global_0, global_5=arima_pre$global_5, global_10=arima_pre$global_10, global_15=arima_pre$global_15, index=year)
        
      # ACF & PACF to determine error structure of residuals; data may follow an ARIMA(p,d,0) model if (i) ACF is exponentially decaying or sinusoidal and (ii) there is a significant spike at lag p in the PACF, but none beyond lag p
        gg_tsdisplay(arima_pre, fn_conc_hg, plot_type = "partial")
        
      # ARIMAs
        # For all: p = 1 (likely- check ACF), d = ndiffs (above), q = 0 because no moving average component
        
      # NO RESERVOIR EFFECTS - run for all rivers
        {
        # No dep
          arima_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + pdq(1,1,0)))
          arima_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + pdq(1,1,0)))
          arima_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + pdq(1,1,0)))
          arima_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + pdq(1,1,0)))
        # Regional dep only
          arima_r00_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ regional_0 + pdq(1,1,0)))
          arima_r05_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ regional_5 + pdq(1,1,0)))
          arima_r10_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ regional_10 + pdq(1,1,0)))
          arima_r15_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ regional_15 + pdq(1,1,0)))
        # 
          arima_r00_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + regional_0 + pdq(1,1,0)))
          arima_r05_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + regional_5 + pdq(1,1,0)))
          arima_r10_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + regional_10 + pdq(1,1,0)))
          arima_r15_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + regional_15 + pdq(1,1,0)))
        # 
          arima_r00_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + regional_0 + pdq(1,1,0)))
          arima_r05_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + regional_5 + pdq(1,1,0)))
          arima_r10_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + regional_10 + pdq(1,1,0)))
          arima_r15_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + regional_15 + pdq(1,1,0)))
        #
          arima_r00_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + regional_0 + pdq(1,1,0)))
          arima_r05_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + regional_5 + pdq(1,1,0)))
          arima_r10_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + regional_10 + pdq(1,1,0)))
          arima_r15_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + regional_15 + pdq(1,1,0)))
        # 
          arima_r00_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + regional_0 + pdq(1,1,0)))
          arima_r05_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + regional_5 + pdq(1,1,0)))
          arima_r10_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + regional_10 + pdq(1,1,0)))
          arima_r15_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + regional_15 + pdq(1,1,0)))
        # Global dep only
          arima_g00_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ global_0 + pdq(1,1,0)))
          arima_g05_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ global_5 + pdq(1,1,0)))
          arima_g10_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ global_10 + pdq(1,1,0)))
          arima_g15_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ global_15 + pdq(1,1,0)))
        # 
          arima_g00_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + global_0 + pdq(1,1,0)))
          arima_g05_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + global_5 + pdq(1,1,0)))
          arima_g10_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + global_10 + pdq(1,1,0)))
          arima_g15_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + tap + global_15 + pdq(1,1,0)))
        # 
          arima_g00_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + global_0 + pdq(1,1,0)))
          arima_g05_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + global_5 + pdq(1,1,0)))
          arima_g10_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + global_10 + pdq(1,1,0)))
          arima_g15_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + maat + global_15 + pdq(1,1,0)))
        #
          arima_g00_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + global_0 + pdq(1,1,0)))
          arima_g05_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + global_5 + pdq(1,1,0)))
          arima_g10_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + global_10 + pdq(1,1,0)))
          arima_g15_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + tap + global_15 + pdq(1,1,0)))
        # 
          arima_g00_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + global_0 + pdq(1,1,0)))
          arima_g05_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + global_5 + pdq(1,1,0)))
          arima_g10_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + global_10 + pdq(1,1,0)))
          arima_g15_05 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + global_15 + pdq(1,1,0)))
            
        # Create empty df
          names <- c("river","AICc","BIC","LBstat","LBpval","model","abbrev")
          arima_df1 <- as.data.frame(matrix(ncol=length(names), nrow=11*4))
          names(arima_df1) <- names
          k=2:7
          j=1
          arima_df1[1:nrow(arima_df1),1] <- riv
        # No Hg dep
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_01))[1],"\\."))[1], "arima_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_02))[1],"\\."))[1], "arima_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_03))[1],"\\."))[1], "arima_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_04))[1],"\\."))[1], "arima_04"); j <- j+1
        # Regional Hg dep
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r00_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r00_01))[1],"\\."))[1], "arima_r00_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r05_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r05_01))[1],"\\."))[1], "arima_r05_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r10_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r10_01))[1],"\\."))[1], "arima_r10_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r15_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r15_01))[1],"\\."))[1], "arima_r15_01"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r00_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r00_02))[1],"\\."))[1], "arima_r00_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r05_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r05_02))[1],"\\."))[1], "arima_r05_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r10_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r10_02))[1],"\\."))[1], "arima_r10_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r15_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r15_02))[1],"\\."))[1], "arima_r15_02"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r00_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r00_03))[1],"\\."))[1], "arima_r00_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r05_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r05_03))[1],"\\."))[1], "arima_r05_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r10_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r10_03))[1],"\\."))[1], "arima_r10_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r15_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r15_03))[1],"\\."))[1], "arima_r15_03"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r00_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r00_04))[1],"\\."))[1], "arima_r00_04"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r05_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r05_04))[1],"\\."))[1], "arima_r05_04"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r10_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r10_04))[1],"\\."))[1], "arima_r10_04"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r15_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r15_04))[1],"\\."))[1], "arima_r15_04"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r00_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r00_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r00_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r00_05))[1],"\\."))[1], "arima_r00_05"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r05_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r05_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r05_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r05_05))[1],"\\."))[1], "arima_r05_05"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r10_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r10_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r10_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r10_05))[1],"\\."))[1], "arima_r10_05"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_r15_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_r15_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_r15_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_r15_05))[1],"\\."))[1], "arima_r15_05"); j <- j+1
        # Global Hg dep
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g00_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g00_01))[1],"\\."))[1], "arima_g00_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g05_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g05_01))[1],"\\."))[1], "arima_g05_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g10_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g10_01))[1],"\\."))[1], "arima_g10_01"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g15_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g15_01))[1],"\\."))[1], "arima_g15_01"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g00_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g00_02))[1],"\\."))[1], "arima_g00_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g05_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g05_02))[1],"\\."))[1], "arima_g05_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g10_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g10_02))[1],"\\."))[1], "arima_g10_02"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g15_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g15_02))[1],"\\."))[1], "arima_g15_02"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g00_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g00_03))[1],"\\."))[1], "arima_g00_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g05_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g05_03))[1],"\\."))[1], "arima_g05_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g10_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g10_03))[1],"\\."))[1], "arima_g10_03"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g15_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g15_03))[1],"\\."))[1], "arima_g15_03"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g00_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g00_04))[1],"\\."))[1], "arima_g00_04"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g05_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g05_04))[1],"\\."))[1], "arima_g05_04"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g10_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g10_04))[1],"\\."))[1], "arima_g10_04"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g15_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g15_04))[1],"\\."))[1], "arima_g15_04"); j <- j+1
        #
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g00_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g00_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g00_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g00_05))[1],"\\."))[1], "arima_g00_05"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g05_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g05_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g05_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g05_05))[1],"\\."))[1], "arima_g05_05"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g10_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g10_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g10_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g10_05))[1],"\\."))[1], "arima_g10_05"); j <- j+1
          arima_df1[j,k] <- c(round(as.numeric(glance(arima_g15_05) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_g15_05) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_g15_05) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_g15_05))[1],"\\."))[1], "arima_g15_05"); j <- j+1
        }
          
      # RESERVOIR EFFECTS - run for Ob' & Yenisey ONLY
        {
        # No dep + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_20 + maat + tap + pdq(1,1,0)))
          arima_res20_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_20 + maat + pdq(1,1,0)))
          arima_res20_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_20 + tap + pdq(1,1,0)))
          arima_res20_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_20 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_r00_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ regional_0 + pdq(1,1,0)))
          arima_res20_r05_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ regional_5 + pdq(1,1,0)))
          arima_res20_r10_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ regional_10 + pdq(1,1,0)))
          arima_res20_r15_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_r00_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + regional_0 + pdq(1,1,0)))
          arima_res20_r05_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + regional_5 + pdq(1,1,0)))
          arima_res20_r10_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + regional_10 + pdq(1,1,0)))
          arima_res20_r15_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_r00_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + regional_0 + pdq(1,1,0)))
          arima_res20_r05_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + regional_5 + pdq(1,1,0)))
          arima_res20_r10_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + regional_10 + pdq(1,1,0)))
          arima_res20_r15_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_r00_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + regional_0 + pdq(1,1,0)))
          arima_res20_r05_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + regional_5 + pdq(1,1,0)))
          arima_res20_r10_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + regional_10 + pdq(1,1,0)))
          arima_res20_r15_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_r00_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + regional_0 + pdq(1,1,0)))
          arima_res20_r05_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + regional_5 + pdq(1,1,0)))
          arima_res20_r10_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + regional_10 + pdq(1,1,0)))
          arima_res20_r15_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + regional_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_g00_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ global_0 + pdq(1,1,0)))
          arima_res20_g05_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ global_5 + pdq(1,1,0)))
          arima_res20_g10_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ global_10 + pdq(1,1,0)))
          arima_res20_g15_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_g00_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + global_0 + pdq(1,1,0)))
          arima_res20_g05_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + global_5 + pdq(1,1,0)))
          arima_res20_g10_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + global_10 + pdq(1,1,0)))
          arima_res20_g15_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + tap + global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_g00_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + global_0 + pdq(1,1,0)))
          arima_res20_g05_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + global_5 + pdq(1,1,0)))
          arima_res20_g10_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + global_10 + pdq(1,1,0)))
          arima_res20_g15_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + maat + global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_g00_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + global_0 + pdq(1,1,0)))
          arima_res20_g05_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + global_5 + pdq(1,1,0)))
          arima_res20_g10_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + global_10 + pdq(1,1,0)))
          arima_res20_g15_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + tap + global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 20 year lag)
          arima_res20_g00_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + global_0 + pdq(1,1,0)))
          arima_res20_g05_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + global_5 + pdq(1,1,0)))
          arima_res20_g10_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + global_10 + pdq(1,1,0)))
          arima_res20_g15_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_20 ~ fn_conc_tss + global_15 + pdq(1,1,0)))
        # No dep + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_01 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_25 + maat + tap + pdq(1,1,0)))
          arima_res25_02 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_25 + maat + pdq(1,1,0)))
          arima_res25_03 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_25 + tap + pdq(1,1,0)))
          arima_res25_04 <- arima_pre %>% model(ARIMA(fn_conc_hg ~ fn_conc_tss + res_cum_cap_25 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_r00_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ regional_0 + pdq(1,1,0)))
          arima_res25_r05_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ regional_5 + pdq(1,1,0)))
          arima_res25_r10_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ regional_10 + pdq(1,1,0)))
          arima_res25_r15_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_r00_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + regional_0 + pdq(1,1,0)))
          arima_res25_r05_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + regional_5 + pdq(1,1,0)))
          arima_res25_r10_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + regional_10 + pdq(1,1,0)))
          arima_res25_r15_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_r00_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + regional_0 + pdq(1,1,0)))
          arima_res25_r05_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + regional_5 + pdq(1,1,0)))
          arima_res25_r10_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + regional_10 + pdq(1,1,0)))
          arima_res25_r15_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_r00_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + regional_0 + pdq(1,1,0)))
          arima_res25_r05_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + regional_5 + pdq(1,1,0)))
          arima_res25_r10_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + regional_10 + pdq(1,1,0)))
          arima_res25_r15_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + regional_15 + pdq(1,1,0)))
        # Regional + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_r00_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + regional_0 + pdq(1,1,0)))
          arima_res25_r05_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + regional_5 + pdq(1,1,0)))
          arima_res25_r10_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + regional_10 + pdq(1,1,0)))
          arima_res25_r15_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + regional_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_g00_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ global_0 + pdq(1,1,0)))
          arima_res25_g05_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ global_5 + pdq(1,1,0)))
          arima_res25_g10_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ global_10 + pdq(1,1,0)))
          arima_res25_g15_01 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_g00_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + global_0 + pdq(1,1,0)))
          arima_res25_g05_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + global_5 + pdq(1,1,0)))
          arima_res25_g10_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + global_10 + pdq(1,1,0)))
          arima_res25_g15_02 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + tap + global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_g00_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + global_0 + pdq(1,1,0)))
          arima_res25_g05_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + global_5 + pdq(1,1,0)))
          arima_res25_g10_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + global_10 + pdq(1,1,0)))
          arima_res25_g15_03 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + maat + global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_g00_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + global_0 + pdq(1,1,0)))
          arima_res25_g05_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + global_5 + pdq(1,1,0)))
          arima_res25_g10_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + global_10 + pdq(1,1,0)))
          arima_res25_g15_04 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + tap + global_15 + pdq(1,1,0)))
        # Global + reservoirs (cumulative capacity, 25 year lag)
          arima_res25_g00_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + global_0 + pdq(1,1,0)))
          arima_res25_g05_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + global_5 + pdq(1,1,0)))
          arima_res25_g10_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + global_10 + pdq(1,1,0)))
          arima_res25_g15_05 <- arima_pre %>% model(ARIMA(fn_conc_hg + res_cum_cap_25 ~ fn_conc_tss + global_15 + pdq(1,1,0)))
          
        # Create empty df
          names <- c("river","AICc","BIC","LBstat","LBpval","model","abbrev")
          arima_df2 <- as.data.frame(matrix(ncol=length(names), nrow=11*4*2))
          names(arima_df2) <- names
          k=2:7
          j=1
          arima_df2[1:nrow(arima_df2),1] <- riv
        # No Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_01))[1],"\\."))[1], "arima_res20_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_02))[1],"\\."))[1], "arima_res20_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_03))[1],"\\."))[1], "arima_res20_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_04))[1],"\\."))[1], "arima_res20_04"); j <- j+1
        # Regional Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r00_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r00_01))[1],"\\."))[1], "arima_res20_r00_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r05_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r05_01))[1],"\\."))[1], "arima_res20_r05_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r10_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r10_01))[1],"\\."))[1], "arima_res20_r10_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r15_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r15_01))[1],"\\."))[1], "arima_res20_r15_01"); j <- j+1
        # Regional Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r00_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r00_02))[1],"\\."))[1], "arima_res20_r00_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r05_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r05_02))[1],"\\."))[1], "arima_res20_r05_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r10_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r10_02))[1],"\\."))[1], "arima_res20_r10_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r15_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r15_02))[1],"\\."))[1], "arima_res20_r15_02"); j <- j+1
        # Regional Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r00_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r00_03))[1],"\\."))[1], "arima_res20_r00_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r05_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r05_03))[1],"\\."))[1], "arima_res20_r05_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r10_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r10_03))[1],"\\."))[1], "arima_res20_r10_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r15_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r15_03))[1],"\\."))[1], "arima_res20_r15_03"); j <- j+1
        # Regional Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r00_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r00_04))[1],"\\."))[1], "arima_res20_r00_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r05_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r05_04))[1],"\\."))[1], "arima_res20_r05_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r10_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r10_04))[1],"\\."))[1], "arima_res20_r10_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_r15_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_r15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_r15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_r15_04))[1],"\\."))[1], "arima_res20_r15_04"); j <- j+1
        # Global Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g00_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g00_01))[1],"\\."))[1], "arima_res20_g00_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g05_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g05_01))[1],"\\."))[1], "arima_res20_g05_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g10_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g10_01))[1],"\\."))[1], "arima_res20_g10_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g15_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g15_01))[1],"\\."))[1], "arima_res20_g15_01"); j <- j+1
        # Global Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g00_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g00_02))[1],"\\."))[1], "arima_res20_g00_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g05_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g05_02))[1],"\\."))[1], "arima_res20_g05_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g10_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g10_02))[1],"\\."))[1], "arima_res20_g10_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g15_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g15_02))[1],"\\."))[1], "arima_res20_g15_02"); j <- j+1
        # Global Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g00_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g00_03))[1],"\\."))[1], "arima_res20_g00_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g05_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g05_03))[1],"\\."))[1], "arima_res20_g05_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g10_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g10_03))[1],"\\."))[1], "arima_res20_g10_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g15_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g15_03))[1],"\\."))[1], "arima_res20_g15_03"); j <- j+1
        # Global Hg dep w/ reservoirs 20 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g00_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g00_04))[1],"\\."))[1], "arima_res20_g00_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g05_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g05_04))[1],"\\."))[1], "arima_res20_g05_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g10_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g10_04))[1],"\\."))[1], "arima_res20_g10_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res20_g15_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res20_g15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res20_g15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res20_g15_04))[1],"\\."))[1], "arima_res20_g15_04"); j <- j+1
        # No Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_01))[1],"\\."))[1], "arima_res25_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_02))[1],"\\."))[1], "arima_res25_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_03))[1],"\\."))[1], "arima_res25_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_04))[1],"\\."))[1], "arima_res25_04"); j <- j+1
        # Regional Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r00_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r00_01))[1],"\\."))[1], "arima_res25_r00_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r05_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r05_01))[1],"\\."))[1], "arima_res25_r05_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r10_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r10_01))[1],"\\."))[1], "arima_res25_r10_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r15_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r15_01))[1],"\\."))[1], "arima_res25_r15_01"); j <- j+1
        # Regional Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r00_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r00_02))[1],"\\."))[1], "arima_res25_r00_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r05_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r05_02))[1],"\\."))[1], "arima_res25_r05_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r10_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r10_02))[1],"\\."))[1], "arima_res25_r10_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r15_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r15_02))[1],"\\."))[1], "arima_res25_r15_02"); j <- j+1
        # Regional Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r00_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r00_03))[1],"\\."))[1], "arima_res25_r00_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r05_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r05_03))[1],"\\."))[1], "arima_res25_r05_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r10_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r10_03))[1],"\\."))[1], "arima_res25_r10_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r15_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r15_03))[1],"\\."))[1], "arima_res25_r15_03"); j <- j+1
        # Regional Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r00_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r00_04))[1],"\\."))[1], "arima_res25_r00_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r05_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r05_04))[1],"\\."))[1], "arima_res25_r05_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r10_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r10_04))[1],"\\."))[1], "arima_res25_r10_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_r15_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_r15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_r15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_r15_04))[1],"\\."))[1], "arima_res25_r15_04"); j <- j+1
        # Global Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g00_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g00_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g00_01))[1],"\\."))[1], "arima_res25_g00_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g05_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g05_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g05_01))[1],"\\."))[1], "arima_res25_g05_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g10_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g10_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g10_01))[1],"\\."))[1], "arima_res25_g10_01"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g15_01) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g15_01) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g15_01))[1],"\\."))[1], "arima_res25_g15_01"); j <- j+1
        # Global Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g00_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g00_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g00_02))[1],"\\."))[1], "arima_res25_g00_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g05_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g05_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g05_02))[1],"\\."))[1], "arima_res25_g05_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g10_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g10_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g10_02))[1],"\\."))[1], "arima_res25_g10_02"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g15_02) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g15_02) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g15_02))[1],"\\."))[1], "arima_res25_g15_02"); j <- j+1
        # Global Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g00_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g00_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g00_03))[1],"\\."))[1], "arima_res25_g00_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g05_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g05_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g05_03))[1],"\\."))[1], "arima_res25_g05_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g10_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g10_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g10_03))[1],"\\."))[1], "arima_res25_g10_03"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g15_03) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g15_03) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g15_03))[1],"\\."))[1], "arima_res25_g15_03"); j <- j+1
        # Global Hg dep w/ reservoirs 25 year lag
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g00_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g00_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g00_04))[1],"\\."))[1], "arima_res25_g00_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g05_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g05_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g05_04))[1],"\\."))[1], "arima_res25_g05_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g10_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g10_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g10_04))[1],"\\."))[1], "arima_res25_g10_04"); j <- j+1
          arima_df2[j,k] <- c(round(as.numeric(glance(arima_res25_g15_04) %>% select(AICc, BIC)),4), round(as.numeric(unlist(augment(arima_res25_g15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[2]),4), round(as.numeric(unlist(augment(arima_res25_g15_04) %>% features(.innov, ljung_box, lag=10, dof=0))[3]),4), unlist(strsplit(names(unlist(arima_res25_g15_04))[1],"\\."))[1], "arima_res25_g15_04"); j <- j+1
        }
        
      # Results
        # Store 
          arima_df1$AICc <- as.numeric(arima_df1$AICc); arima_df1$BIC <- as.numeric(arima_df1$BIC); arima_df1$LBstat <- as.numeric(arima_df1$LBstat); arima_df1$LBpval <- as.numeric(arima_df1$LBpval)
          arima_df2$AICc <- as.numeric(arima_df2$AICc); arima_df2$BIC <- as.numeric(arima_df2$BIC); arima_df2$LBstat <- as.numeric(arima_df2$LBstat); arima_df2$LBpval <- as.numeric(arima_df2$LBpval)
        # Summarize dataframes
          dim(arima_df1); dim(na.omit(arima_df1))
          dim(arima_df2); dim(na.omit(arima_df2))
        # Store df1, for all rivers except Ob' & Yenisey
          arima_df <- arima_df1
        # Bind both df1 and df2 (for Ob' & Yenisey only)
          arima_df <- rbind(arima_df1, arima_df2)
        # Order df by increasing BIC
          arima_df <- arima_df[order(arima_df$BIC, decreasing=F),]
        # Summary of final df
          dim(arima_df); dim(na.omit(arima_df))
          head(arima_df); summary(arima_df)
        # Show top models
          arima_df[1:3,] # Top 3 models, based on BIC
          arima_df[1:10,] # Top 10
          #arima_df %>% slice(which.min(BIC)) %>% select(BIC, model) # lowest BIC; for models with low number of observations, AIC can be biased and thus AICc is preferred (https://otexts.com/fpp3/selecting-predictors.html
        
      # Final candidate model: QA/QC
        # Optionally export results
          #write.csv(arima_df, paste0(dir,riv,"_DRM_results.csv"))
        # Store candidate model
          arima_final <- arima_02
        # Plot residuals. Should resemble white noise- uncorrelated, mean of zero and ideally constant variance (homoscedastic) and normally distributed) | https://otexts.com/fpp3/diagnostics.html
          arima_final %>% gg_tsresiduals()
            # Raw data for residuals by year, if desired
              #as_tibble(residuals(arima_final, type = "regression")) # Nt
              #as_tibble(residuals(arima_final, type = "innovation")) # Et
        # Ljung-Box test- p > 0.05 indicates autocorrelations are indistinguishable from a series of white noise
          augment(arima_final) %>% features(.innov, ljung_box, lag=10, dof=0) # Use lag=10 for non-seasonal data, dof (i.e., K) = 0 if autocorrelations calculated from raw data (rather than residuals) (Hyndman & Athanasopoulos, 2021)
        # Summary  
          report(arima_final)
  
      # Model summary
        # Coefficients: for covariates of standard regression (Yt) and error structure(s) (Nt)
        # s.e. = std error of coeffs
        # Nt = regression residuals
        # Et = ARIMA residuals: sigma^2 = std deviation of white noise, i.e. Et, where √sigma^2 = Et; normally & independently distributed (NID) around 0
        # Yt = intercept + coeffs*predictors + coeffs*Nt + Et
        # log likelihood: logarithm of the probability of the observed data coming from the estimated model. Given (p,d,q), loglik is maximized by ARIMA when finding parameter estimates
        # For writing model w/ coefficients, see https://otexts.com/fpp3/non-seasonal-arima.html & https://otexts.com/fpp3/regarima.html


  # OTHER ####
    
    # Clean global environment
      {rm(arima_01,
      arima_02,
      arima_03,
      arima_04,
      arima_r00_01,
      arima_r05_01,
      arima_r10_01,
      arima_r15_01,
      arima_r00_02,
      arima_r05_02,
      arima_r10_02,
      arima_r15_02,
      arima_r00_03,
      arima_r05_03,
      arima_r10_03,
      arima_r15_03,
      arima_r00_04,
      arima_r05_04,
      arima_r10_04,
      arima_r15_04,
      arima_r00_05,
      arima_r05_05,
      arima_r10_05,
      arima_r15_05,
      arima_g00_01,
      arima_g05_01,
      arima_g10_01,
      arima_g15_01,
      arima_g00_02,
      arima_g05_02,
      arima_g10_02,
      arima_g15_02,
      arima_g00_03,
      arima_g05_03,
      arima_g10_03,
      arima_g15_03,
      arima_g00_04,
      arima_g05_04,
      arima_g10_04,
      arima_g15_04,
      arima_g00_05,
      arima_g05_05,
      arima_g10_05,
      arima_g15_05,
      arima_res20_01,
      arima_res20_02,
      arima_res20_03,
      arima_res20_04,
      arima_res20_r00_01,
      arima_res20_r05_01,
      arima_res20_r10_01,
      arima_res20_r15_01,
      arima_res20_r00_02,
      arima_res20_r05_02,
      arima_res20_r10_02,
      arima_res20_r15_02,
      arima_res20_r00_03,
      arima_res20_r05_03,
      arima_res20_r10_03,
      arima_res20_r15_03,
      arima_res20_r00_04,
      arima_res20_r05_04,
      arima_res20_r10_04,
      arima_res20_r15_04,
      arima_res20_r00_05,
      arima_res20_r05_05,
      arima_res20_r10_05,
      arima_res20_r15_05,
      arima_res20_g00_01,
      arima_res20_g05_01,
      arima_res20_g10_01,
      arima_res20_g15_01,
      arima_res20_g00_02,
      arima_res20_g05_02,
      arima_res20_g10_02,
      arima_res20_g15_02,
      arima_res20_g00_03,
      arima_res20_g05_03,
      arima_res20_g10_03,
      arima_res20_g15_03,
      arima_res20_g00_04,
      arima_res20_g05_04,
      arima_res20_g10_04,
      arima_res20_g15_04,
      arima_res20_g00_05,
      arima_res20_g05_05,
      arima_res20_g10_05,
      arima_res20_g15_05,
      arima_res25_01,
      arima_res25_02,
      arima_res25_03,
      arima_res25_04,
      arima_res25_r00_01,
      arima_res25_r05_01,
      arima_res25_r10_01,
      arima_res25_r15_01,
      arima_res25_r00_02,
      arima_res25_r05_02,
      arima_res25_r10_02,
      arima_res25_r15_02,
      arima_res25_r00_03,
      arima_res25_r05_03,
      arima_res25_r10_03,
      arima_res25_r15_03,
      arima_res25_r00_04,
      arima_res25_r05_04,
      arima_res25_r10_04,
      arima_res25_r15_04,
      arima_res25_r00_05,
      arima_res25_r05_05,
      arima_res25_r10_05,
      arima_res25_r15_05,
      arima_res25_g00_01,
      arima_res25_g05_01,
      arima_res25_g10_01,
      arima_res25_g15_01,
      arima_res25_g00_02,
      arima_res25_g05_02,
      arima_res25_g10_02,
      arima_res25_g15_02,
      arima_res25_g00_03,
      arima_res25_g05_03,
      arima_res25_g10_03,
      arima_res25_g15_03,
      arima_res25_g00_04,
      arima_res25_g05_04,
      arima_res25_g10_04,
      arima_res25_g15_04,
      arima_res25_g00_05,
      arima_res25_g05_05,
      arima_res25_g10_05,
      arima_res25_g15_05)}
        