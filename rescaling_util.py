from xarray import open_dataset, Dataset, DataArray, where, ones_like
from utils import latname, lonname, precipname, tempname
import numpy as np
from xesmf import Regridder
from os.path import isfile
from decimal import *
from square_mask import square_mask

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

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

def mult_scaling(month, forcing, target_name, clim_loc, lats, lons, latmask, lonmask, adjust, Unique_ID):

    if lats.size == 1:
        lats, latmask = np.array([lats]), np.array(latmask)
    if lons.size == 1:
        lons, lonmask = np.array([lons]), np.array(lonmask) 
        
    rounded_lats = np.array([int(10**5*i) for i in lats])/10**5
    rounded_lons = np.array([int(10**5*i) for i in lats])/10**5

    reference_DA = DataArray(np.ones((lats.size, lons.size), dtype=np.float64),
                             dims = ['latitude', 'longitude'],
                             coords = {'latitude':rounded_lats,
                                       'longitude':rounded_lons
                                      })
    
    if (target_name is None) | (adjust == 'neither'):
        scaled_t2m = ones_like(reference_DA)[latmask, lonmask].values
        scaled_tp = ones_like(reference_DA)[latmask, lonmask].values
        
    else:
        
    
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

        # lims
        latmin = np.min(lats)
        latmax = np.max(lats)

        # limits
        clim = clim.where((clim.latitude >= latmin) & (clim.latitude <= latmax), drop=True)
        target = target.where((target.latitude >= latmin) & (target.latitude <= latmax), drop=True)

        # get bounds 
        gb = grid_bounds(clim)

        weights_loc = '/users/jk/20/achereque/SnowProjects2/src/B-TIM/weights/'
        weights_name = forcing+'r_target_'+target_name+'_'+str(month).zfill(2)+'.nc'

        reuse_weights = False

        if (isfile(weights_loc+'bil_'+weights_name)) & (reuse_weights):
            regridder_bil = Regridder(grid_bounds(target), grid_bounds(clim), 
                                      'bilinear', periodic=True, 
                                      weights=weights_loc+'bil_'+weights_name, reuse_weights=True)
        else:
            regridder_bil = Regridder(grid_bounds(target), grid_bounds(clim), 
                                      'bilinear', periodic=True)

            regridder_bil.to_netcdf(weights_loc+'bil_'+weights_name)

        if (isfile(weights_loc+'cons_'+weights_name)) & (reuse_weights):
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

        out = out.where(out.latitude.isin(lats), drop=True)
        clim = clim.where(clim.latitude.isin(lats), drop=True)

        # order
        latdiff = lats[1] - lats[0]
        londiff = lons[1] - lons[0]

        if latdiff > 0:
            out = out.sortby('latitude')
            clim = clim.sortby('latitude')
        else:
            out = out.sortby('latitude', ascending=False)
            clim = clim.sortby('latitude', ascending=False)

        if londiff > 0:
            out = out.sortby('longitude')
            clim = clim.sortby('longitude')
        else:
            out = out.sortby('longitude', ascending=False)
            clim = clim.sortby('longitude', ascending=False)


        print(reference_DA.longitude.values, 
              clim.t2m.longitude.values, 
              target.t2m.longitude.values,
              out.t2m.longitude.values)

        if (adjust=='both') | (adjust=='t2m'):
            fraction_t2m = where((out.t2m/clim.t2m) >= 99/100, out.t2m/clim.t2m, 99/100)
            fraction_t2m = where((fraction_t2m) <= 100/99, fraction_t2m, 100/99)
            fraction_t2m = fraction_t2m.fillna(1.)

            scaled_t2m =  (reference_DA * fraction_t2m)[latmask,lonmask].values
        else:
            scaled_t2m = ones_like(reference_DA)[latmask, lonmask].values

        if (adjust=='both') | (adjust=='tp'):

            fraction_tp = where((out.tp/clim.tp) >= 1/4, out.tp/clim.tp, 1/4)
            fraction_tp = where(fraction_tp <= 4, fraction_tp, 4)

            fraction_tp = where((clim.tp < 0.01)|(out.tp < 0.01), 1., fraction_tp)

            scaled_tp =  (reference_DA * fraction_tp)[latmask, lonmask].values
        else:
            scaled_tp = ones_like(reference_DA)[latmask, lonmask].values

        if adjust == 'both':
            fig, axs = map_ax(ncols=2)

            levs = axs[0].contourf(lons[lonmask], lats[latmask], 
                                   scaled_t2m, 
                                   transform=ccrs.PlateCarree(), vmin=0.985, 
                                   vmax=1.02, levels=28, extend='neither', cmap='RdBu_r')
            plt.colorbar(levs, ax=axs[0])

            levs2 = axs[1].contourf(lons[lonmask], lats[latmask], 
                                    scaled_tp, transform=ccrs.PlateCarree(), 
                                    vmin=0, vmax=2, extend='neither', levels=28, cmap='RdBu_r')
            plt.colorbar(levs2, ax=axs[1])

            plt.savefig('/users/jk/20/achereque/SnowProjects2/src/B-TIM/snapshot_scaling_factors/'+Unique_ID+'_'+str(month).zfill(2)+'.png')
            plt.close()
    
    return scaled_t2m, scaled_tp
