from numpy import isin
from os import listdir
from netCDF4 import Dataset
from xarray import open_mfdataset

month_names_aug = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'March', 'April', 'May', 'June', 'July']

t2m_freq = {'JRA55': 8, 'ERA5':24, 'MERRA2':24}
tp_freq = {'JRA55': 8, 'ERA5':24, 'MERRA2':24}
latname = {'JRA55':'latitude', 'ERA5':'latitude', 'MERRA2':'lat'}
lonname = {'JRA55':'longitude', 'ERA5':'longitude', 'MERRA2':'lon'}
precipname = {'JRA55':'tp', 'ERA5':'tp', 'MERRA2':'PRECTOTLAND'}
tempname = {'JRA55':'t2m', 'ERA5':'t2m', 'MERRA2':'T2M'}


def len_month(m, year, leapday=True):
    '''Takes in month (Jan = 1) and year, returns the number of days in the month.'''
    if isin(m, [4, 6, 9, 11]):
        return 30
    elif m == 2:
        if (leapday) & (year % 100 == 0) & (year % 400 != 0):
            return 28
        elif (leapday) & (year % 4 == 0):
            return 29
        else:
            return 28
    else:
        return 31
    

def prepare_filenames(forcing, data_loc, snow_season):
    
    prefix = lambda var,m,y: data_loc+var+'/'+forcing+'/'+forcing+'_'+var+'_'+str(m).zfill(2)+'_'+str(y)+'/'
    
    tp_AugDec = [[prefix('tp',i,snow_season[0])+n for n in listdir(prefix('tp',i,snow_season[0])) if '.nc' in n] for i in range(8, 13)]
    tp_JanJul = [[prefix('tp',i,snow_season[1])+n for n in listdir(prefix('tp',i,snow_season[1])) if '.nc' in n] for i in range(1,8)]

    t2m_AugDec = [[prefix('t2m',i,snow_season[0])+n for n in listdir(prefix('t2m',i,snow_season[0])) if '.nc' in n] for i in range(8, 13)]
    t2m_JanJul = [[prefix('t2m',i,snow_season[1])+n for n in listdir(prefix('t2m',i,snow_season[1])) if '.nc' in n] for i in range(1,8)]
    
    tp_files = tp_AugDec + tp_JanJul
    t2m_files = t2m_AugDec + t2m_JanJul
    
    return t2m_files, tp_files

def read_month(forcing, t2m_fname, tp_fname, month, y):
    
    if forcing == 'MERRA2':
        
        t2m = open_mfdataset(t2m_fname, combine='by_coords')
        tp = open_mfdataset(tp_fname, combine='by_coords')

        t2m = t2m.rename({latname[forcing]:'latitude',
                          lonname[forcing]:'longitude'})
        tp = tp.rename({latname[forcing]:'latitude',
                        lonname[forcing]:'longitude'})
        
        full_lat, full_lon = t2m['latitude'].values, t2m['longitude'].values
 
    else: 
        if len(t2m_fname) == 1:
            t2m = Dataset(t2m_fname[0])
        if len(tp_fname) == 1:
            tp = Dataset(tp_fname[0])

        full_lat, full_lon = t2m[latname[forcing]][:], t2m[lonname[forcing]][:]
    
    return full_lat, full_lon, t2m, tp

def read_day(forcing, data, forcing_var, step, latmask, lonmask):
    
    decode_var = {'tp':precipname[forcing], 't2m':tempname[forcing]}
    
    if forcing == 'MERRA2':
        output = data[decode_var[forcing_var]][step,latmask,lonmask].values
    else:
        output = data[decode_var[forcing_var]][step,latmask,lonmask]

    return output

def standardize_precip(forcing, pr_freq, t2m_freq, data):
    '''Convert into [m per precipitation time step] from native units.'''
    
    if forcing == 'JRA55':
        data = data / (pr_freq * 1000) #[m/tp step from mm/day]
        data = (pr_freq/t2m_freq) * data #[m/t2m step] 
    elif forcing == 'MERRA2':
        data = 3600 * data / 1000 #[m/precip step] from [mm/s]
        data = (pr_freq/t2m_freq) * data #[m/t2m step]
    elif forcing == 'ERA5':
        data = (pr_freq/t2m_freq) * data #[m/t2m step]
   
    return data  
    
def monthly_out_name(Unique_ID, month, mixed_pr, year_tag):
    
    savename = Unique_ID + '_forced_swe_' + month + '_'

    if mixed_pr:
        savename += 'mixedpr_'
        
    savename = savename + year_tag + '.nc'
    
    return savename

