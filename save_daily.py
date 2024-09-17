from xarray import Dataset, DataArray
from numpy import moveaxis

def save_daily(lats, lons, times, snf_record, density_record, out_fname):
    
    if lats.size == 1:
        lats = lats[0]
    if lons.size == 1:
        lons = lons[0]
    
    sdepDA = DataArray(data = moveaxis(snf_record, (0,1,2), (1,2,0)),
                        dims = ['time','lat','lon'],
                        coords = {
                            'time': times,
                            'lat': (['lat'], lats),
                            'lon':(['lon'], lons)
                            
                        },
                        attrs = {
                            'description': 'snow depth in metres of snow',
                            'units': 'm',
                            'standard_name': 'surface_snow_thickness'
                        }
                      )

    sdenDA = DataArray(data = moveaxis(density_record, (0,1,2), (1,2,0)),
                       dims = ['time','lat','lon'],
                       coords = {
                            'time': times,
                            'lat': (['lat'], lats),
                            'lon':(['lon'], lons),
                       },
                       attrs = {
                            'description': 'snow density in kilograms per cubic metre',
                            'units': 'kg/m3',
                            'standard_name': 'surface_snow_density'
                       }
                      )
                       
    
    dataset = Dataset({'snow_depth': sdepDA,
                       'density': sdenDA}
                     )

    dataset.lon.attrs = {'units':'degrees_east', 'long_name':'longitude'}
    dataset.lat.attrs = {'units':'degrees_north', 'long_name':'latitude'}
    
    dataset.to_netcdf(out_fname)
    dataset.close()
    print('saved to netcdf:', out_fname)
