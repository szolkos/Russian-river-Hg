// This code and associated Earth Engine assets can be accessed at: https://code.earthengine.google.com/?scriptPath=users%2Fszolkos%2FRussianRivHg%3Aera5

// Script written by G. Fiske, Woodwell Climate Research Center, adapted by S. Zolkos

// Add in the data
var year1 = 1979;
var year2 = 2020;
var years = ee.List.sequence(year1, year2);
var era5_2mt = ee.ImageCollection('ECMWF/ERA5/MONTHLY')
            //.select('mean_2m_air_temperature') // for air temp
            .select('total_precipitation') // for precip
            .filterDate(year1 + '-01-01',year2 + '-12-31')
            .map(function(img)\{
              var d = ee.Date(ee.Number(img.get('system:time_start')));
              var y = ee.Number(d.get('year'));
              //img = img.subtract(273.15); // convert K to C; turn off for precip
              return img.set(\{'year':y\});
            });
var era5_2mt_Count = era5_2mt.size();
print('collection count: ', era5_2mt_Count);

// aggregate by year
var byYear = ee.ImageCollection.fromImages(
      years.map(function(y)\{
            return era5_2mt.filterMetadata('year', 'equals', y)
                        //.select('mean_2m_air_temperature')
                        //.mean().round() // <== takes the mean of all the images per year, for temp
                        //.rename('meanTemp') // <== for temp
                        .select('total_precipitation')
                        .sum()  // <== sum all images per year, for precip
                        .rename('sumPrecip') // <== for precip
                        .set('year', y)
                        .set('date', ee.Date.fromYMD(y,1,1).millis());
      })
    );
print('yearly count', byYear.size());
print('byYear', byYear);

// add to map
Map.centerObject(rusriv, 3.8);
var palette = ['#EFEFEF', '#F3D4C9', '#F7B7A2', '#F89980','#F07B77','#e75b6e','#CC4B60','#AD4052','#903544','#732A36'];
Map.addLayer(byYear.sort('year', false).first().clip(pechora),\{min:-15,max:15, palette: palette\},"Mean (most recent year)");
Map.addLayer(rusriv.style(\{color: 'black', width: 1, fillColor: '00000000'\}), \{\}, 'Russian river watersheds');

// make a chart of the mean temperature values by year
// For precip change 'meanTemp' to 'sumPrecip'
// use mean() or stdDev()
var pChart = ui.Chart.image.seriesByRegion(byYear, pechora, ee.Reducer.mean(), 'sumPrecip', 1000, 'year', 'name') // change meanTemp to sumPrecip if need be
      .setChartType('ScatterChart')
      .setOptions(\{
        //title: 'mean yearly temperature',
        title: 'total annual precipitation',
        lineWidth: 1,
        pointSize: 3
      });
print(pChart);}
