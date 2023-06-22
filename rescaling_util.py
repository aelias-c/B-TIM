from xarray import open_dataset, Dataset, DataArray, where, ones_like
from utils import latname, lonname, precipname, tempname
import numpy as np
from xesmf import Regridder
from os.path import isfile
from square_mask import square_mask

def grid_bounds(data, pwr_10 = 5):

    fac = 10**pwr_10
    latitudes = np.array([int(fac * lat) for lat in data['latitude']])
    longitudes = np.array([int(fac * lon) for lon in data['longitude']])
    
    lonstep = (longitudes[1] - longitudes[0])
    latstep = (latitudes[1] - latitudes[0])
    
    lat_b = np.arange(latitudes[0]-latstep, latitudes[-1]+latstep, latstep)
    lon_b = np.arange(longitudes[0]-lonstep, longitudes[-1]+lonstep, lonstep)

    lat = latitudes - (latstep)/2
    lon = longitudes - (lonstep)/2

    return {'lon': lon/fac, 'lat': lat/fac, 'lon_b': lon_b/fac, 'lat_b': lat_b/fac}

def translate_0_360(data, pwr_10 = 5):
    
    fac = 10**pwr_10
    lons = np.array([int(fac * l) for l in data['longitude']])
    lons_0_360 = lons % (fac * 360)
    
    data = data.assign_coords(longitude = lons_0_360 / fac)
    
    return data.sortby('longitude').sortby('latitude')
    
def translate_180(data, pwr_10 = 5):
        
    fac = 10**pwr_10
    lons = np.array([int(fac * l) + (fac * 180) for l in data['longitude']])
    lons_180 = (lons % (fac * 360)) - (fac * 180)
    
    data = data.assign_coords(longitude = lons_180 / fac)
    
    return data.sortby('longitude').sortby('latitude')

def mult_scaling(month, forcing, target_name, clim_loc, lats, lons, latmask, lonmask, adjust, Unique_ID, weights_loc = './weights/'):
    
    clim = open_dataset(clim_loc+forcing+'_'+str(month).zfill(2)+'_mm.nc')
    target = open_dataset(clim_loc+target_name+'_'+str(month).zfill(2)+'_mm.nc')

    original_lons = np.sum(clim.longitude<0)
    if original_lons != 0: #if some longitudes are negative
        clim = translate_180(clim)
        target = translate_180(target)
    else:
        clim = translate_0_360(clim)
        target = translate_0_360(target)

    # transpose
    clim = clim.transpose('latitude', 'longitude')
    target = target.transpose('latitude', 'longitude')
    
    # get bounds 
    gb = grid_bounds(clim)
    
    weights_name = forcing+'r_target_'+target_name+'_'+str(month).zfill(2)+'.nc'
    
    if isfile(weights_loc+'bil_'+weights_name):
        regridder_bil = Regridder(grid_bounds(target), grid_bounds(clim), 
                                  'bilinear', periodic=True, 
                                  weights=weights_loc+'bil_'+weights_name, reuse_weights=True)
    else:
        regridder_bil = Regridder(grid_bounds(target), grid_bounds(clim), 
                                  'bilinear', periodic=True)
        
        regridder_bil.to_netcdf(weights_loc+'bil_'+weights_name)
        
    if isfile(weights_loc+'cons_'+weights_name):
        regridder_cons = Regridder(grid_bounds(target), grid_bounds(clim), 
                                   'conservative', weights=weights_loc+'cons_'+weights_name, 
                                   reuse_weights=True)
    else:
        regridder_cons = Regridder(grid_bounds(target), grid_bounds(clim), 
                                   'conservative')
        
        regridder_cons.to_netcdf(weights_loc+'cons_'+weights_name)

    # regrid data
    out = Dataset({'tp':DataArray(data = regridder_cons(target.tp.values), 
                                      dims = ['latitude', 'longitude'], 
                                      coords = {'latitude': gb['lat_b'][1:],
                                                'longitude': gb['lon_b'][1:]}),
                   't2m':DataArray(data = regridder_bil(target.t2m.values), 
                                      dims = ['latitude', 'longitude'], 
                                      coords = {'latitude': gb['lat_b'][1:],
                                                'longitude': gb['lon_b'][1:]})})

    out = out.where(out.latitude.isin(lats))
    clim = clim.where(clim.latitude.isin(lats))
    
    if (adjust=='both') | (adjust=='t2m'):
        fraction_t2m = where(out.t2m/clim.t2m>= 9/10, out.t2m/clim.t2m, 9/10)
        fraction_t2m = where(fraction_t2m <= 10/9, fraction_t2m, 10/9)
        fraction_t2m = fraction_t2m.fillna(1.)
        
        scaled_t2m =  fraction_t2m[latmask,lonmask].values
    else:
        scaled_t2m = ones_like(clim.t2m[latmask,lonmask]).values
        
    if (adjust=='both') | (adjust=='tp'):
        clim_tp = where((clim.tp < 1e-3)|(out.tp < 1e-3), np.nan, clim.tp)
        out_tp = where((clim.tp < 1e-3)|(out.tp < 1e-3), np.nan, out.tp)
        
        fraction_tp = where(out_tp/clim_tp >= 1/4, out_tp/clim_tp, 1/4)
        fraction_tp = where(fraction_tp <= 4, fraction_tp, 4)
        fraction_tp = fraction_tp.fillna(1.)
        
        scaled_tp =  fraction_tp[latmask, lonmask].values
    else:
        scaled_tp = ones_like(clim.tp[latmask,lonmask]).values
        
        
    return scaled_t2m, scaled_tp
