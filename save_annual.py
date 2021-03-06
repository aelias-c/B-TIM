from xarray import Dataset, DataArray

def save_annual(experiment_name, output_loc, mixed_pr, year_tag, lats, lons, 
                ptot_record, sftot_record, SWEmax_record):
    
    ptotDA = DataArray(data = ptot_record,
                    dims = ['latitude', 'longitude'],
                    coords = {
                        'latitude': (['latitude'], lats),
                        'longitude':(['longitude'], lons)
                    },
                    attrs = {
                        'description': 'total precipitation (frozen and liquid)',
                        'units': "mm",
                        'standard_name': 'lwe_thickness_of_precipitation_amount'
                    }
                  )
    sftotDA = DataArray(data = sftot_record,
                    dims = ['latitude', 'longitude'],
                    coords = {
                        'latitude': (['latitude'], lats),
                        'longitude':(['longitude'], lons)
                    },
                    attrs = {
                        'description': 'total snowfall (sum of lwe falling when 2m-temperature is below freezing)',
                        'units': "mm",
                        'standard_name': 'lwe_thickness_of_surface_snow_amount'
                    }
                  )
    SWEmaxDA = DataArray(data = SWEmax_record,
                    dims = ['latitude', 'longitude'],
                    coords = {
                        'latitude': (['latitude'], lats),
                        'longitude':(['longitude'], lons)
                    },
                    attrs = {
                        'description': 'hydrological-year maximum snow water equivalent',
                        'units': "mm"
                    }
                  )
    
    output_dataset = Dataset(
        {
            'ptot': ptotDA,
            'sftot': sftotDA, 
            'swemax': SWEmaxDA 
        }
    )

    savename = experiment_name + '_forced_'
    if mixed_pr:
        savename += 'mixedpr_'
    savename += year_tag + '.nc'

    output_dataset.to_netcdf(output_loc + savename)
    
    output_dataset.close()