import numpy as np
import pandas as pd
import xarray as xr
import rasterio
import matplotlib.pyplot as plt
import xesmf as xe


def make_catchments_dataset():
    '''
    
    '''
    # --------------------------------------------
    # Part 0. open tif with catchment masks
    # --------------------------------------------
    data_path = './Data/'
    catchments = rasterio.open(data_path+'riv.tif').read(1)
    fine = xr.open_mfdataset(data_path+'IF.nc')

    # --------------------------------------------
    # Part 1. convert catchments to xarray Dataset
    # --------------------------------------------
    lon = np.linspace(35.5,171.5,1281)
    lat = np.linspace(45,73,251)
    catchment_names = {1: 'Onega', 
                       2: 'Mezen', 
                       3: 'Pechora', 
                       4: 'North Dvina', 
                       5: 'Kolyma', 
                       6: 'Lena', 
                       7: 'Ob', 
                       8: 'Yenisey', 
                       255: None}
    
    catchments_all = xr.DataArray(catchments, 
                                  coords={'lat':(lat[0:-1]+lat[1:])/2, 
                                          'lon':(lon[0:-1]+lon[1:])/2}, 
                                  dims=['lat','lon'])

    catchments = xr.Dataset()
    for i in catchment_names:
        name = catchment_names[i]
        # replace nans with 0 in .where() function
        # this is important to conserving area in downscaling
        catchments[name] = catchments_all.where(catchments_all==i, 0)/i
    
    catchments = catchments.drop(None)
    
    # --------------------------------------------
    # Part 2. regrid catchment masks to 2 x 2.5 to match IFs
    # --------------------------------------------
    lat_res = fine.lat
    lon_res = fine.lon
    ds_out = xr.Dataset({'lat':(['lat'], lat_res), 
                         'lon':(['lon'], lon_res)})
    regridder = xe.Regridder(catchments, ds_out, 'bilinear', )#periodic=True)
    regridder  # print basic regridder information.
    catchments_original = catchments.copy() # make copy for visual comparison below
    # perform regridding
    catchments = regridder(catchments)
    catchments = catchments.assign_coords({'area':fine['area']})

    return catchments


def read_streets():
    '''
    
    '''
    data_path = './Data/'
    tmp = pd.read_csv(data_path+'500_year_ALW.csv', header=1)
    streets_regions = ['NA','SA','EU','USSR','AF','AS','OC','Tot']
    categories = ['Air', 'Air_Hg0', 'Air_Hg2', 'Land_Water', 'Total_Hg']

    column_headers = []
    for i in streets_regions:
        for j in categories:
            column_headers.append(j+'_'+i)
    tmp.columns = ['Year']+column_headers

    yrs = np.linspace(1510, 2010, 501)
    df = pd.DataFrame({'Year': yrs})
    df = df.merge(tmp, left_on=['Year'], right_on=['Year'], how='outer')
    df = df.interpolate(method='linear', limit_direction='forward')

    return df

def compile_emissions():
    '''
    
    '''
    geo_Hg0 = 250.
    ocean_Hg0 = 3836.
    soil_Hg0 = 853.
    bb_Hg0 = 321.
    
    
    # read streets 500 years
    data_path = './Data/'
    tmp = pd.read_csv(data_path+'500_year_ALW.csv', header=1)
    streets_regions = ['NA','SA','EU','USSR','AF','AS','OC','Tot']
    categories = ['Air', 'Air_Hg0', 'Air_Hg2', 'Land_Water', 'Total_Hg']

    column_headers = []
    for i in streets_regions:
        for j in categories:
            column_headers.append(j+'_'+i)
    tmp.columns = ['Year']+column_headers

    yrs = np.linspace(1510, 2010, 501)
    streets = pd.DataFrame({'Year': yrs})
    streets = streets.merge(tmp, left_on=['Year'], right_on=['Year'], how='outer')
    streets = streets.interpolate(method='linear', limit_direction='forward')
    
    anthro_emissions = pd.DataFrame({
        'Year':streets.Year, 
        'AF_Hg0':(streets.Air_Hg0_AF+streets.Air_Hg0_OC),
        'AF_HgII':(streets.Air_Hg2_AF+streets.Air_Hg2_OC),
        'AS_Hg0':(streets.Air_Hg0_AS),
        'AS_HgII':(streets.Air_Hg2_AS),                        
        'EU_Hg0':(streets.Air_Hg0_EU+streets.Air_Hg0_USSR),
        'EU_HgII':(streets.Air_Hg2_EU+streets.Air_Hg2_USSR),
        'NA_Hg0':streets.Air_Hg0_NA,
        'NA_HgII':streets.Air_Hg2_NA,
        'SA_Hg0':streets.Air_Hg0_SA,
        'SA_HgII':streets.Air_Hg2_SA,
    })
    
    non_anthro_emissions = pd.DataFrame({
        'Year':streets.Year, 
        'Geogenic_Hg0': np.repeat(geo_Hg0,len(streets.Year)),
        'Ocean_Hg0':    np.repeat(ocean_Hg0,len(streets.Year)),
        'Soil_Hg0':     np.repeat(soil_Hg0,len(streets.Year)),
        'Biomass_Burning_Hg0': np.repeat(bb_Hg0,len(streets.Year)),})
    
    # subset just last 50 years of record
    all_emissions = anthro_emissions.merge(non_anthro_emissions, how='outer', on='Year').iloc[-51:]
    
    ## mutliply Hg0:HgII speciation from 2010 by Streets et al. (2019) 2010-2015 trends
    Streets_addition = pd.read_csv(data_path+'Anthro_2010-2015.csv')
    Streets_regions_dict = {
        'AF':['Northern Africa', 'Western Africa', 'Eastern Africa', 'Southern Africa', 'Middle East', 'Oceania'],
        'AS':['South Asia', 'East Asia', 'Southeast Asia', 'Japan'],
        'EU':['OECD Europe', 'Eastern Europe', 'Former USSR'], 
        'NA':['Canada','USA','Central America'],
        'SA':['South America']}

    update = pd.DataFrame()
    for i in Streets_regions_dict.keys():
        update[i] = np.array([0,0,0,0,0,0])
        for j in Streets_regions_dict[i]:

            update[i] += Streets_addition[Streets_addition['World region']== j].iloc[0,2:].values
    update = update.T.copy()
    update.columns = ([2010,2011,2012,2013,2014,2015])

    df = pd.DataFrame({'Year':[2010.,2011.,2012.,2013.,2014.,2015.]})
    for i in ['AF','AS','EU','NA','SA']:
        Hg0 = all_emissions[f'{i}_Hg0'].iloc[-1]
        Hg2 = all_emissions[f'{i}_HgII'].iloc[-1]
        f_Hg0 = Hg0/(Hg0+Hg2)

        tmp_Hg0 = f_Hg0*update.loc[i].values
        tmp_HgII = (1-f_Hg0)*update.loc[i].values

        df[f'{i}_Hg0'] = tmp_Hg0
        df[f'{i}_HgII'] = tmp_HgII

    df['Geogenic_Hg0'] = np.repeat(geo_Hg0,6)
    df['Ocean_Hg0'] = np.repeat(ocean_Hg0,6)
    df['Soil_Hg0'] = np.repeat(soil_Hg0,6)
    df['Biomass_Burning_Hg0'] = np.repeat(bb_Hg0,6)

    all_emissions = all_emissions.append(df.iloc[1:,:], ignore_index=True)

    for i in all_emissions.columns:
        all_emissions[i] = all_emissions[i].astype(float)

    return all_emissions

def make_deposition_fields(all_emissions, fine, fine_speciated):
    '''
    
    '''
    emission_magnitudes = xr.Dataset()
    for column in all_emissions.columns[1:]:
        emission_magnitudes[column] = xr.DataArray(data   = all_emissions[column],
                                                   dims   = ['time'],
                                                   coords = {'time': all_emissions['Year']})
    emission_magnitudes = emission_magnitudes.sel(time=slice(1960,2015))

    ## Make total deposition Dataset
    deposition = xr.Dataset()
    init_var = 'AF_Hg0'
    deposition['EU_dep'] = fine[init_var]*emission_magnitudes[init_var]*0
    for i in ['EU_Hg0','EU_HgII']:
        deposition['EU_dep'] += fine[i]*emission_magnitudes[i]

    deposition['natural_dep'] = fine[init_var]*emission_magnitudes[init_var]*0
    for i in ['Geogenic_Hg0', 'Ocean_Hg0','Soil_Hg0','Biomass_Burning_Hg0']:
        deposition['natural_dep'] += fine[i]*emission_magnitudes[i]

    deposition['all_dep'] = fine[init_var]*emission_magnitudes[init_var]*0
    for i in emission_magnitudes.data_vars:
        deposition['all_dep'] += fine[i]*emission_magnitudes[i]

    deposition_monthly = deposition.copy()
    deposition = deposition.sum('month')

    ## Make speciated deposition Dataset (equivalent of above)
    deposition_speciated = xr.Dataset()
    init_var = 'AF_Hg0'
    for j in ['Hg2wet','Hg2dry','Hg0dry','Hg0uptake']:
        deposition_speciated[f'{j}_all_dep'] = fine_speciated[f'{j}_{init_var}']*emission_magnitudes[init_var]*0

    for i in emission_magnitudes.data_vars:
        for j in ['Hg2wet','Hg2dry','Hg0dry','Hg0uptake']:
            deposition_speciated[f'{j}_all_dep'] += fine_speciated[f'{j}_{i}']*emission_magnitudes[i]

    deposition_speciated_monthly = deposition_speciated.copy()
    deposition_speciated = deposition_speciated.sum('month')
    
    return deposition, deposition_speciated

def make_deposition_table(deposition, catchments):
    '''
    
    '''
    deposition_table = pd.DataFrame({
        'River': ['Onega','North Dvina', 'Mezen','Pechora','Ob','Yenisey','Lena','Kolyma'],
        'Start Year':[1980, 1980, 1980, 1983, 1979, 1980, 2004, 1982],
        'End Year':  [1992, 2001, 2002, 2002, 2008, 2011, 2011, 1987]})

    deltas_abs = []
    deltas_pct = []
    intercepts = []
    dep_2015 = []
    dep_1970 = []
    time_series = pd.DataFrame({'Year':deposition.time})

    for i in deposition_table['River']:

        start = deposition_table[deposition_table['River'] == i]['Start Year'].values.item()
        end = deposition_table[deposition_table['River'] == i]['End Year'].values.item()
        river = i

        time_series[i] = (deposition['all_dep']*catchments[i]).sum(('lat','lon')).values
        modern = (deposition['all_dep']*catchments[i]).sum(('lat','lon')).sel(time=2015.).values
        historical = (deposition['all_dep']*catchments[i]).sum(('lat','lon')).sel(time=1970.).values
        start_dep = (deposition['all_dep']*catchments[i]).sum(('lat','lon')).sel(time=start).values
        end_dep = (deposition['all_dep']*catchments[i]).sum(('lat','lon')).sel(time=end).values
        delta = end_dep - start_dep
        delta_pct = ( (end_dep-start_dep)/start_dep )*1e2

        deltas_abs.append(delta)
        deltas_pct.append(delta_pct)
        dep_2015.append(modern)
        dep_1970.append(historical)
        intercepts.append(start_dep)

    deposition_table['Intercept'] = intercepts
    deposition_table['Intercept'] = deposition_table['Intercept'].astype(float)
    deposition_table['$\Delta$ (total change)'] = deltas_abs
    deposition_table['$\Delta$ (%)'] = deltas_pct
    deposition_table['Dep. 1970 (Mg yr-1)'] = dep_1970
    deposition_table['Dep. 2015 (Mg yr-1)'] = dep_2015
    deposition_table['Dep. 1970 (Mg yr-1)'] = deposition_table['Dep. 1970 (Mg yr-1)'].astype(float)
    deposition_table['Dep. 2015 (Mg yr-1)'] = deposition_table['Dep. 2015 (Mg yr-1)'].astype(float)
    pd.options.display.float_format = "{:,.1f}".format
    
    return deposition_table

def get_timeseries_speciated(catchments, deposition_speciated):
    '''
    
    '''
    rivers = ['Onega','North Dvina','Mezen','Pechora','Ob','Yenisey','Lena','Kolyma']
    time_series_speciated = pd.DataFrame({'Year':deposition.time})
    for river in rivers:
        for dep_type in deposition_speciated.data_vars:
            time_series_speciated[f'{river}_{dep_type}'] = (deposition_speciated[dep_type]*catchments[river]).sum(('lat','lon')).values
    return time_series_speciated

def get_timeseries_all(catchments, deposition):
    '''
    
    '''
    time_series_all = pd.DataFrame({'Year':deposition.time})
    for i in catchments.data_vars:
        river = i
        time_series_all[i] = (deposition['all_dep']*catchments[i]).sum(('lat','lon')).values
    return time_series_all

def get_timeseries_EU(catchments, deposition):
    '''
    
    '''
    time_series_EU = pd.DataFrame({'Year':deposition.time})
    for i in catchments.data_vars:
        river = i
        time_series_EU[i] = (deposition['EU_dep']*catchments[i]).sum(('lat','lon')).values
    return time_series_EU

