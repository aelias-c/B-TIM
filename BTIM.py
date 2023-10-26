import numpy as np
from pandas import date_range

from utils import len_month, monthly_out_name, month_names_aug, prepare_filenames, read_month, read_day, t2m_freq, tp_freq, standardize_precip, standardize_temp
from square_mask import square_mask
from time_step import Brasnett
from save_daily import save_daily
from save_annual import save_annual

import CONFIG as cfg

# ------- Initialize -------- #

snow_season = [cfg.year, cfg.year+1]
t2m_files, tp_files = prepare_filenames(cfg.forcing, cfg.data_loc, snow_season)

year_tag = str(snow_season[0])+'_'+str(snow_season[1])
print('Processing snow year: ' + year_tag)

# --- Step month by month --- #  
for i,m in enumerate([8,9,10,11,12,1,2,3,4,5,6,7]):
    
    current_y = snow_season[0]
    if m < 8:
        current_y = snow_season[1]
        
    days_in_month = len_month(m, current_y, cfg.leapdays)
    
    full_lat, full_lon, t2m, pr = read_month(m, cfg.forcing, t2m_files[i], tp_files[i])
                
    if i == 0:
        
        latmask = square_mask(full_lat, full_lon, cfg.latminmax)
        lonmask = np.ones_like(full_lon, dtype='bool')
        
        nlats, nlons = np.sum(latmask), np.sum(lonmask)
        lats, lons = full_lat[latmask], full_lon[lonmask]
        
        if lats.size == 1:
            lats, latmask = np.array([lats]), np.array([latmask])
        if lons.size == 1:
            lons, lonmask = np.array([lons]), np.array([lonmask])


        # ------ Set up records for the year once ------ #
        ptot_record = np.zeros((nlats, nlons)) #[m], total precip 
        sftot_record = np.zeros((nlats, nlons)) #[m water equivalent], total snowfall
        SWEmax_record = np.zeros((nlats, nlons)) #[mm water equivalent], maximum SWE

        # - Set up prognostic variable grids only once - #
        old_depth = np.zeros((nlats, nlons)) #[m snow depth]
        old_dens = np.zeros((nlats, nlons)) #[kg/m3]
   
    t2m_scale, tp_scale = np.ones((nlats, nlons)), np.ones((nlats, nlons))
    
    # ----- Set up daily records for the month ----- #
    snf_record = np.zeros((nlats, nlons, days_in_month)) #[m snow], snow depth
    density_record = np.zeros((nlats, nlons, days_in_month)) #[kg/m3], snow density
    
    # ------------- Step through month ------------- #
    day = 0
    for step in range(days_in_month * t2m_freq[cfg.forcing]):
        
         # ------------ Read in precip data ------------- #
        t2m_steps_per_pr = t2m_freq[cfg.forcing] // tp_freq[cfg.forcing]
        prate = read_day(cfg.forcing, pr, 'tp', step // t2m_steps_per_pr,  latmask, lonmask)
        prate = tp_scale * standardize_precip(cfg.forcing, 
                                              tp_freq[cfg.forcing], 
                                              t2m_freq[cfg.forcing], 
                                              prate) #[m water] per precipitation time step
        prate[prate < 0] = 0
        ptot_record += prate  
        
        # ---------- Read in temperature data ---------- #
        if (step == 0) & (i == 0):
            # initially use the same values for t2m_last as for t2m_air
            read_t2m = read_day(cfg.forcing, t2m, 't2m', step, latmask, lonmask)
            t2m_last = t2m_scale * standardize_temp(cfg.forcing, read_t2m) #[K]
        else:
            t2m_last = t2m_air #[K]
        read_t2m = read_day(cfg.forcing, t2m, 't2m', step, latmask, lonmask)
        t2m_air = t2m_scale * standardize_temp(cfg.forcing, read_t2m) #[K]
        TSFC = np.stack((t2m_last, t2m_air)) - 273.15 #[degrees C]  
        tavg = np.mean(TSFC, axis=0) #[degrees C]
        
        # ------ Record snowfall where tavg < 0C ------- #
        sftot_record[tavg <= 0] += prate[tavg <= 0]
        
        # --------- Linearly interpolate temp ---------- #
        hours_per_step = 24 // t2m_freq[cfg.forcing]
        T_hr = np.ones((hours_per_step+1, TSFC.shape[1], TSFC.shape[2]))
        for hr in range(hours_per_step+1):
            T_hr[hr,:,:] = (TSFC[1,:,:] - TSFC[0,:,:]) * hr / hours_per_step + TSFC[0,:,:] #[degrees C]
        
        # -------- Calculate mean hourly precip -------- #
        TP_hr = prate / hours_per_step #[m] in one hour
        
        # ----------- Time-step by one chunk ----------- #
        old_depth, old_dens, swe = Brasnett(cfg.mixed_pr, T_hr, TP_hr, old_depth, old_dens) #[m], [kg/m3], [mm]

        # --------- Track any record-high SWE ---------- #
        SWEmax_record = np.maximum(SWEmax_record, swe) #[mm water equivalent]

        # --- Record daily depth and density values ---- #
        if (step + 1) % t2m_freq[cfg.forcing] == 0: #last time step each day
            print('day ', day)
            snf_record[:,:,day] = old_depth #[m snow]
            density_record[:,:,day] = old_dens #[kg/m3]
            
            # --- Record daily depth and density values ---- #
            
            # ----------- Write to monthly file ------------ #
            if (day + 1) == days_in_month: #last time step of last day of month
                times = date_range(str(current_y)+'-'+str(m).zfill(2)+'-'+'01', periods=days_in_month, freq='D')
                
                #set up save name according to settings
                out_fname = cfg.output_loc + monthly_out_name(cfg.Unique_ID, 
                                                                month_names_aug[i], 
                                                                cfg.mixed_pr,
                                                                year_tag)
                save_daily(lats, lons, times, 
                           snf_record, density_record, 
                           out_fname)
                
            day += 1

    if i != 11:
        pr.close()
        t2m.close()      
        
# ------ Save accumulated records to file ------ #
save_annual(cfg.Unique_ID, cfg.output_loc, cfg.mixed_pr, 
            year_tag, lats, lons, 
            ptot_record, sftot_record, SWEmax_record)

pr.close()
t2m.close()
