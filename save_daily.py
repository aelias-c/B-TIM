from xarray import Dataset, DataArray

def save_daily(lats, lons, times, snf_record, density_record, out_fname):
    
    sdepDA = DataArray(data = snf_record,
                        dims = ['time', 'latitude', 'longitude'],
                        coords = {
                            'latitude': (['latitude'], lats),
                            'longitude':(['longitude'], lons),
                            'time': times
                        },
                        attrs = {
                            'description': 'snow depth in centimetres of snow',
                            'units': "cm",
                            'standard_name': 'surface_snow_thickness'
                        }
                      )

    sdenDA = DataArray(data = density_record,
                    dims = ['time', 'latitude', 'longitude'],
                    coords = {
                        'latitude': (['latitude'], lats),
                        'longitude':(['longitude'], lons),
                        'time': times
                    },
                    attrs = {
                        'description': 'snow density in kilograms per cubic metre',
                        'units': 'kg/m3',
                        'standard_name': 'surface_snow_density'
                    }
                  )
                       
    
    dataset = Dataset({'snow_depth': sdepDA,
                       'density': sdenDA
                      }
                     )

    dataset.to_netcdf(out_fname)
    dataset.close()
    print('saved to netcdf:', out_fname)
