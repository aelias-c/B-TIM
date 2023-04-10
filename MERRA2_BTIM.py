from netCDF4 import Dataset
from os import listdir
from xarray import open_mfdataset

forcing = 'MERRA2'
data_loc = '/users/jk/20/achereque/TEST_DELETE/B-TIM/forcing/'
t2m_freq = 24
pr_freq = 24
temp_in_K = True
pr_in_mhr = False
lon_0_360 = False

tp_name = 'PRECTOTLAND'
t2m_name = 'T2M'
latitude_name = 'lat'
longitude_name = 'lon'


def prepare_filenames(snow_season):
    
    prefix = lambda var,m,y: data_loc+forcing+'_'+var+'_'+str(m).zfill(2)+'_'+str(y)+'/'
    
    
    tp_AugDec = [[prefix('tp',i,snow_season[0])+n for n in listdir(prefix('tp',i,snow_season[0])) if '.nc' in n] for i in range(8, 13)]
    tp_JanJul = [[prefix('tp',i,snow_season[1])+n for n in listdir(prefix('tp',i,snow_season[1])) if '.nc' in n] for i in range(1,8)]

    t2m_AugDec = [[prefix('t2m',i,snow_season[0])+n for n in listdir(prefix('t2m',i,snow_season[0])) if '.nc' in n] for i in range(8, 13)]
    t2m_JanJul = [[prefix('t2m',i,snow_season[1])+n for n in listdir(prefix('t2m',i,snow_season[1])) if '.nc' in n] for i in range(1,8)]
    
    tp_files = tp_AugDec + tp_JanJul
    t2m_files = t2m_AugDec + t2m_JanJul
    
    return t2m_files, tp_files

def read_month(t2m_file, tp_file, month, y):
    t2m = open_mfdataset(t2m_file, combine='by_coords')
    tp = open_mfdataset(tp_file, combine='by_coords')
    
    t2m = t2m.rename({latitude_name:'latitude',
                      longitude_name:'longitude'})
    tp = tp.rename({latitude_name:'latitude',
                    longitude_name:'longitude'})

    full_lat, full_lon = t2m['latitude'].values, t2m['longitude'].values
    
    return full_lat, full_lon, t2m, tp

def read_day(data, var, step, latmask, lonmask):
    decode_var = {'tp':tp_name, 't2m':t2m_name}
    output = data[decode_var[var]][step,latmask,lonmask].values

    return output
        

def standardize_precip(data):
    '''Convert into m per precipitation time step from native units'''
    data = 3600 * data / 1000 #[m/precip step] from [mm/s]
    return data  
