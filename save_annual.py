from xarray import Dataset, DataArray

def save_annual(Unique_ID, output_loc, mixed_pr, year_tag, lats, lons, 
                ptot_record, sftot_record, SWEmax_record):
    
    if lats.size == 1:
        lats = lats[0]
    if lons.size == 1:
        lons = lons[0]
    
    ptotDA = DataArray(data = ptot_record,
                    dims = ['latitude', 'longitude'],
                    coords = {
                        'latitude': (['latitude'], lats),
                        'longitude':(['longitude'], lons)
                    },
                    attrs = {
                        'description': 'total precipitation (frozen and liquid)',
                        'units': "m",
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
                        'units': "m",
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
                        'description': 'water year maximum snow water equivalent',
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

    savename = Unique_ID + '_forced_'
    if mixed_pr[0] != mixed_pr[1]:
        savename += 'mixedpr_'
    savename += year_tag + '.nc'

    output_dataset.to_netcdf(output_loc + savename)
    
    output_dataset.close()
