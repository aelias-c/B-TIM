from netCDF4 import Dataset

forcing = 'ERA5'
data_loc = 'forcing/'
chunks = 24
pr_freq = 1
temp_in_K = True
pr_in_mhr = True
lon_0_360 = True

tp_name = 'tp'
t2m_name = 't2m'
latitude_name = 'latitude'
longitude_name = 'longitude'


def prepare_filenames(snow_season):
    prefix = data_loc + forcing
    tp_AugDec = [prefix + '_tp_' + str(f'{i:02}') + '_' + str(snow_season[0]) + '.nc' for i in range(8, 13)]
    tp_JanJul = [prefix + '_tp_'+str(f'{i:02}') + '_' + str(snow_season[1]) + '.nc' for i in range(1, 8)]

    t2m_AugDec = [prefix + '_t2m_' + str(f'{i:02}') + '_' + str(snow_season[0]) + '.nc' for i in range(8, 13)]
    t2m_JanJul = [prefix + '_t2m_' + str(f'{i:02}') + '_' + str(snow_season[1]) + '.nc' for i in range(1, 8)]
    
    tp_files = tp_AugDec + tp_JanJul
    t2m_files = t2m_AugDec + t2m_JanJul
    
    return t2m_files, tp_files

def read_month(t2m_file, tp_file, month, y):
    t2m = Dataset(t2m_file)
    tp = Dataset(tp_file)
    
    full_lat, full_lon = t2m[latitude_name][:], t2m[longitude_name][:]
    
    return full_lat, full_lon, t2m, tp

def read_day(data, var, step, latmask, lonmask):
    if var == 'tp':
        v = tp_name
    elif var == 't2m':
        v = t2m_name
        
    return data[v][step,latmask,lonmask]

def standardize_precip(data):
    '''Convert into m per precipitation time step from native units'''
    return data


        
