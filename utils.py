from numpy import isin
from os import listdir
from netCDF4 import Dataset
from xarray import open_mfdataset

month_names_aug = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'March', 
                   'April', 'May', 'June', 'July']

#number of time steps per day in forcing data
t2m_freq = {'JRA55': 8, 'ERA5':24, 'MERRA2':24} 
tp_freq = {'JRA55': 8, 'ERA5':24, 'MERRA2':24}

# if unsure of names of variables or coordinates, use ncdump
latname = {'JRA55':'latitude', 'ERA5':'latitude', 'MERRA2':'lat'}
lonname = {'JRA55':'longitude', 'ERA5':'longitude', 'MERRA2':'lon'}

precipname = {'JRA55':'tp', 'ERA5':'tp', 'MERRA2':'PRECTOTLAND'}
tempname = {'JRA55':'t2m', 'ERA5':'t2m', 'MERRA2':'T2M'}

def len_month(m, year, leapday=True):
    '''Takes in month and year, returns the number of days in the month.
    
    Args:
        m (int): month number (Jan = 1)
        year (int): year, used in combination with month
        leapday (bool): whether to account for leapdays in February
        
    Returns:
        len_month (int): number of days in m,year
    '''
    
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
    '''Gathers all temperature and precipitation filenames for one snow season.
    
    Args:
        forcing (str): name of forcing dataset
        data_loc (str): path to directory with forcing data
        snow_season (tuple): list of the form [YYYY, YYYY+1] containing the 
            calendar years covered by the snow season being modeled.
        
    Returns:
        t2m_files (list): 12 filenames or lists of filenames with temperature 
            forcing data, beginning with August.
        tp_files (list): 12 filenames or lists of filenames with precipitation 
            forcing data, beginning with August.
    '''
    
    tp_AugDec = [data_loc+forcing+'_tp_'+str(m).zfill(2)+'_'+str(snow_season[0])+'.nc' for m in range(8,13)]
    tp_JanJul = [data_loc+forcing+'_tp_'+str(m).zfill(2)+'_'+str(snow_season[1])+'.nc' for m in range(1,8)]
    
    t2m_AugDec = [data_loc+forcing+'_t2m_'+str(m).zfill(2)+'_'+str(snow_season[0])+'.nc' for m in range(8,13)]
    t2m_JanJul = [data_loc+forcing+'_t2m_'+str(m).zfill(2)+'_'+str(snow_season[1])+'.nc' for m in range(1,8)]

    tp_files = tp_AugDec + tp_JanJul
    t2m_files = t2m_AugDec + t2m_JanJul
    
    return t2m_files, tp_files

def read_month(month, forcing, t2m_fname, tp_fname):
    '''Extracts coordinate and variable data for one month into xarray dataset 
       or netCDF4 dataset. 
    
    Args:
        month (int): month to read, Jan = 1
        forcing (str): name of forcing dataset
        t2m_fname (str or list of str): file name(s) containing temperature
            data for the month.
        tp_fname (str or list of str): file name(s) containing precipitation
            data for the month.
        
    Returns:
        full_lat (ndarray): latitudes extracted from file
        full_lon (ndarray): longitudes extracted from file
        t2m (dataset): temperature data for month
        tp (dataset): precipitation data for month
    '''
    
    #add forcings here when there are multiple .nc files per month
    if isin(forcing, ['MERRA2']): 
        
        t2m = open_mfdataset(t2m_fname, combine='by_coords')
        tp = open_mfdataset(tp_fname, combine='by_coords')

        t2m = t2m.rename({latname[forcing]:'latitude',
                          lonname[forcing]:'longitude'})
        tp = tp.rename({latname[forcing]:'latitude',
                        lonname[forcing]:'longitude'})
        
        full_lat, full_lon = t2m['latitude'].values, t2m['longitude'].values
 
    else: #add forcings here when there is only one .nc file per month
        if (len(tp_fname) == 1)&(len(t2m_fname) == 1):
            t2m = Dataset(t2m_fname[0])
            tp = Dataset(tp_fname[0])
        else:
            t2m = Dataset(t2m_fname)
            tp = Dataset(tp_fname)

        full_lat, full_lon = t2m[latname[forcing]][:], t2m[lonname[forcing]][:]
    
    return full_lat, full_lon, t2m, tp

def read_day(forcing, data, forcing_var, step, latmask, lonmask):
    '''Extracts forcing data for just one time step and region.
    
    Args:
        forcing (str): name of forcing dataset
        data (dataset)
        forcing_var (str): name of variable
        step (int): time step in month to extract
        latmask (ndarray): boolean mask for region
        lonmask (ndarray): boolean mask for region
        
    Returns:
        output (ndarray)
    '''
    
    decode_var = {'tp':precipname[forcing], 't2m':tempname[forcing]}
    
    if isin(forcing, ['MERRA2']):
        
        if (latmask.size == 1) & (lonmask.size == 1):
            output = data[decode_var[forcing_var]][step].values
        elif (latmask.size == 1):
            output = data[decode_var[forcing_var]][step,lonmask].values
        elif (lonmask.size == 1):
            output = data[decode_var[forcing_var]][step,latmask].values
        else:
            output = data[decode_var[forcing_var]][step,latmask,lonmask].values
            
    else:
        if (latmask.size == 1) & (lonmask.size == 1):
            output = data[decode_var[forcing_var]][step]
        elif (latmask.size == 1):
            output = data[decode_var[forcing_var]][step,lonmask]
        elif (lonmask.size == 1):
            output = data[decode_var[forcing_var]][step,latmask]
        else:
            output = data[decode_var[forcing_var]][step,latmask,lonmask]

    return output

def standardize_precip(forcing, pr_freq, t2m_freq, data):
    '''Convert into [m per temp time step] from native units.'''
    
    if forcing == 'JRA55':
        data = data / (pr_freq * 1000) #[m/tp step from mm/day]
        data = (pr_freq/t2m_freq) * data #[m/t2m step] 
        
    elif forcing == 'MERRA2':
        data = 3600 * data / 1000 #[m/precip step] from [mm/s]
        data = (pr_freq/t2m_freq) * data #[m/t2m step]
        
    elif forcing == 'ERA5':
        data = (pr_freq/t2m_freq) * data #[m/t2m step]
   
    return data  

def standardize_temp(forcing, data):
    '''Convert into [K] from native units.'''
   
    return data  
    
def monthly_out_name(Unique_ID, month, mixed_pr, year_tag):
    '''Construct output filename.'''
    
    savename = Unique_ID + '_forced_swe_' + month + '_'
   
    if mixed_pr[0] != mixed_pr[1]:
        savename += 'mixedpr_'
        
    savename = savename + year_tag + '.nc'
    
    return savename

