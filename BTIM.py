import sys
import numpy as np
from pandas import date_range

from utils import len_month, monthly_out_name, month_names_aug
from square_mask import square_mask
from time_step import Brasnett
from save_daily import save_daily
from save_annual import save_annual

import config as exd

year = int(sys.argv[1])

snow_season = [year, year+1]
year_tag = str(year)+'_'+str(year+1)

t2m_files, tp_files = exd.prepare_filenames(snow_season)

print('Processing snow year: ' + year_tag)

for i,m in enumerate([8,9,10,11,12,1,2,3,4,5,6,7]):
    
    days_in_month = len_month(m, snow_season[1], exd.leapdays)
    
    full_lat, full_lon, t2m, tp = exd.read_month(t2m_files[i], tp_files[i], m)
    
    if i == 0:
        
        latmask, lonmask = square_mask(full_lat, full_lon, 
                                       exd.latminmax, exd.lonminmax, 
                                       exd.lon_0_360)
        
        nlats, nlons = np.sum(latmask), np.sum(lonmask)
        
        ### Set up yearly records
        ptot_record = np.zeros((nlats, nlons)) #[m], total precip 
        sftot_record = np.zeros((nlats, nlons)) #[m water equivalent], total snowfall
        SWEmax_record = np.zeros((nlats, nlons)) #[m water equivalent], record swe
        
        ### Set up daily records
        snf_record = np.zeros((nlats, nlons, days_in_month))
        density_record = np.zeros((nlats, nlons, days_in_month))
        
        ### Set up arrays for prognostic variables
        old_depth = np.zeros((nlats, nlons)) #[cm snow]
        old_dens = np.zeros((nlats, nlons)) #[kg/m^3]
        
    else:
        
        ### Set up daily records
        snf_record = np.zeros((nlats, nlons, days_in_month))
        density_record = np.zeros((nlats, nlons, days_in_month))
        
    day = 1
    for step in range(days_in_month * exd.t2m_freq):
        if (step == 0) & (i == 0):
            #for first step only, use the same values for 
            #t2m last as for t2m_air
            t2m_last = exd.read_day(t2m, 't2m', 0, latmask, lonmask) #[K]
        else:
            t2m_last = t2m_air #[K]
            
        t2m_air = exd.read_day(t2m, 't2m', 0, latmask, lonmask) #[K]
        tsfc = np.stack((t2m_last, t2m_air)) - 273.15 #[degrees C]
        
        pcpn = exd.read_day(tp, 'tp', step//exd.pr_freq, latmask, lonmask) 
        pcpn = exd.standardize_precip(pcpn) #[mm per precipitation time step]
        pcpn[pcpn < 0] = 0
        ptot_record += pcpn
        
        tavg = np.mean(tsfc, axis=0)
        sftot_record[tavg <= 0] += pcpn[tavg <= 0]
        
        # set up forcing to be hourly
        hours_per_chunk = int(24/exd.pr_freq)
        T_hr = np.ones((hours_per_chunk+1, tsfc.shape[1], tsfc.shape[2]))
        for hr in range(hours_per_chunk+1): #linearly interpolate temperature to be hourly
            T_hr[hr,:,:] = (tsfc[1,:,:] - tsfc[0,:,:]) * hr / hours_per_chunk + tsfc[0,:,:]
        TP_hr = pcpn / hours_per_chunk #[m] in one hour  
    
        # step forward one time step, sub-routine runs hourly
        old_depth, old_dens, swe = Brasnett(hours_per_chunk, exd.mixed_pr, T_hr, TP_hr, old_depth, old_dens)
        
        SWEmax_record = np.maximum(SWEmax_record, swe)
        
        if (step) % exd.t2m_freq == 0:
            print('day: ', day)
            snf_record[:,:,day-1] = old_depth
            density_record[:,:,day-1] = old_dens
            
            if day == days_in_month:
                yr = snow_season[0]
                if m < 8:
                    yr = snow_season[1]
                    
                times = date_range(str(yr)+'-'+str(m).zfill(2)+'-'+'01', periods=days_in_month, freq='D')
                out_fname = exd.output_loc + monthly_out_name(exd.unique_ID, 
                                                              month_names_aug[i], 
                                                              exd.mixed_pr,
                                                              year_tag)
                save_daily(full_lat[latmask], full_lon[lonmask], 
                           times, snf_record, density_record, out_fname)
            
            day += 1
    if i < 11:
        tp.close()
        t2m.close()

save_annual(exd.unique_ID, exd.output_loc, exd.mixed_pr, year_tag, full_lat[latmask], full_lon[lonmask], 
            ptot_record, sftot_record, SWEmax_record)

tp.close()
t2m.close()
        
