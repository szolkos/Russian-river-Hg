# Russian river sediment and mercury
## Introduction
Source code for the manuscript: Multi-decadal declines in particulate mercury and sediment export from Russian rivers in the pan-Arctic basin. 

## Authors
- [Scott Zolkos](https://www.researchgate.net/profile/Scott-Zolkos)
- Alexander V. Zhulidov
- Tatiana Yu. Gurtovaya
- Vyacheslav V. Gordeev
- [Sergey Berdnikov](https://scholar.google.com/citations?user=3zmlVS8AAAAJ&hl=en&oi=sra)
- Nadezhda Pavlova
- Evgenia A. Kalko
- Yana A. Kuklina
- Danil A. Zhulidov
- Lyudmila S. Kosmenko
- [Alexander I. Shiklomanov](https://scholar.google.com/citations?user=7mevZAYAAAAJ&hl=en&oi=sra)
- Anya Suslova
- [Benjamin M. Geyman](https://scholar.google.com/citations?user=qs1jNdQAAAAJ&hl=en&oi=ao)
- [Colin P. Thackray](https://scholar.google.com/citations?user=2Bj2PrIAAAAJ&hl=en&oi=sra)
- [Elsie M. Sunderland](http://bgc.seas.harvard.edu/index.html)
- [Suzanne E. Tank](https://www.suzannetank.net)
- [James W. McClelland](https://scholar.google.com/citations?user=jRd24DQAAAAJ&hl=en)
- [Robert G.M. Spencer](https://www.spencerbiogeochem.org)
- [David P. Krabbenhoft](https://www.researchgate.net/profile/David-Krabbenhoft)
- Richard Robarts
- [Robert M. Holmes](https://arcticgreatrivers.org)

## Data
### The following data used in this manuscript can be downloaded from the Arctic Data Center:
- *concentrations.csv*: Concentrations of total suspended solids (TSS; mg/L) and mercury in TSS (i.e., particulate mercury, PHg; µg/g). Used in subsequent analyses.  
- *Q_Onega.csv*: Onega River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_NorthDvina.csv*: Northern Dvina River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_Mezen.csv*: Mezen River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_Pechora.csv*: Pechora River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_Ob.csv*: Ob' River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_Yenisey.csv*: Yenisey River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_Lena.csv*: Lena River daily discharge data (m^3/s). Used in subsequent analyses.  
- *Q_Kolyma.csv*: Kolyma River daily discharge data (m^3/s). Used in subsequent analyses.  
- *GRanD_Russia.csv*: Data on Russian reservoirs compiled from the GRanD database.  
- *res_cum_cap_w_lags.csv*: Reservoir cumulative capacity (km^3) derived from the GRanD database.  
- *GEOS_Chem_Hg_dep.csv*: Annual atmospheric Hg deposition (mg/y, 1960-2015) to each watershed.  
- *era5_temp_precip_1979_2020.csv*: Mean annual air temperature (MAAT) and total annual precipitation (1979-2020) for each watershed.  
- *lit_comparison.csv*: PHg values from this study and from the literature.  
- *WRTDS_EGRET*: Results from WRTDS analyses of TSS and PHg trends; used in plotting.  
 
## Scripts
### Modeling
*WRTDS.R*: Flow-weighted TSS and PHg concentration and fluxes, and their trends, for each river.  
*GEOS_Chem.R*: Atmospheric Hg deposition magnitudes to each watershed.  

### Data processing
#### *Note: These scripts process the data summarized above, which was compiled as detailed in the Methods of this manuscript. Users of these scripts will first need to update the working directory pathway within script #1.*  
*1_read_data.R*: Load packages, annotated bibliography, parameter definitions, read data (total suspended solids [TSS], particulate mercury [PHg], discharge).  
*2_prep_conc.R*: Prepare concentration data: (1) Calculate mean of replicates; (2) Correct for C-Q spatiotemporal offsets; (3) Calculate PHg.  
*3_wrtds.R*: Weighted regressions on time, discharge, and season (WRTDS) modeling of concentration and flux.  
*4_Hg_dep.R*: Read in data for Hg deposition from GEOS-Chem.  
*5_reservoirs.R*: Data from GRanD (Lehnher et al., 2011) on reservoirs in Russian river watersheds.  
*6a_climate.js*: Google Earth Engine script for deriving mean annual air temperature (MAAT) and total annual precipitation (TAP) by watershed from ERA5 climate reanalysis data.  
*6b_climate.R*: Trends in mean annual air temperature (MAAT) and total annual precipitation (TAP) by watershed.  
*7_drm.R*: Trends in sediment and mercury concentrations and fluxes, and drivers.  
*8_phg_lit_compare.r*:  Compare measured PHg (µg/g) from this study w/ observations in global rivers.  

## Packages
- plotrix
- ggplot2
- gdata
- lubridate
- gridExtra
- ggpubr
- fable
- fpp3
- feasts
- tsfeatures
- magrittr
- GGally
- tibble
- plyr
- dplyr
- tidyr
- ggtext
