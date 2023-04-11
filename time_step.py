import numpy as np

def warm_snow_aging(DENSITY, SWE, timestep, wet_settle_mask):
    
    sdep = (SWE[wet_settle_mask] / DENSITY[wet_settle_mask])  #[cm]
    denmax = 700. - ((204.70 / sdep) * (1 - np.exp(-sdep / 0.673)))
    den_diff = denmax - DENSITY[wet_settle_mask]

    TIMFAC = np.exp(np.log(den_diff[den_diff > 0.1] / 200.) - (2.778e-6 * timestep))

    del_den_warm = np.zeros_like(sdep)
    del_den_warm[den_diff > 0.1] = den_diff[den_diff > 0.1] - 200 * TIMFAC
    
    return del_den_warm

def cold_snow_aging(DENSITY, SWE, TDD, icl, cold_settle_mask):
    C2 = np.where((icl == 2) | (icl == 6), -28/1000, -21/1000)[cold_settle_mask]
    del_den_cold = 2 * np.exp(0.08 * TDD[cold_settle_mask]) * (0.6 * SWE[cold_settle_mask]) * np.exp(C2 * DENSITY[cold_settle_mask]) #[kg/m3]
    
    return del_den_cold     
                      

def temp_melt(TDD, GAMMA, temp_melt_mask):
    '''
    Update depth when temperature exceeds T_melt.
    
    Args:
        TDD (np.array): temperature for sub-step.
        GAMMA (np.array)
        temp_melt_mask (np.array)
    '''
    
    del_SWE = np.zeros_like(TDD)
    del_SWE[temp_melt_mask] = GAMMA[temp_melt_mask] * TDD[temp_melt_mask] #[m water]
    
    return del_SWE

def rain_melt(RAIN, T_diff, rain_melt_mask):    
    '''
    Melt snow wherever there is rain according to rain_melt_mask.
    
    Args:
        RAIN (np.array): total liquid precipitation (m water)
        T_diff (np.array): Temerature for time step in degrees C
        rain_melt_mask
    '''
    
    ### define some constants
    rhow = 1000 #[kg/m^3], density of water
    Cw = 4.18e3 #[J], specific heat of water
    rhoice = 917 #[kg/m^3], density of ice
    Lf = 0.334e6 #[J/kg], latent heat of fusion of water

    ### compute rain 
    ### RAIN has units [m water] equiv. [1 m3 water per m2]
    mass_water_per_m2 = rhow * RAIN[rain_melt_mask] #[kg/m2] 
    heat_from_rain = mass_water_per_m2 * Cw * T_diff[rain_melt_mask] #[J/m2]
    
    del_SWE = np.zeros_like(RAIN)
    del_SWE[rain_melt_mask] = heat_from_rain / (Lf * rhoice) #[m water]

    return del_SWE


def new_snow_density(T):
    '''
    New snow density follows Hedstrom and Pomeroy (1998).
    
    Args:
        T (np.array): Temerature for time step in degrees C
    '''
    
    ### density is temperature dependent
    rhosfall_cold = 67.9 + 51.3 * np.exp(T / 2.6)
    
    ### rhosfall_warm will only be used if mixed precip is active
    rhosfall_warm =  np.minimum(119.2 + (20 * T), 200.)
  
    rhosfall = np.where(T <= 0, rhosfall_cold, rhosfall_warm)
    return rhosfall

def hourly_melt_rate(DENSITY, boreal_mask):
    ''' 
    Find hourly melt rate from daily rate based on 
    Kuusito (1980) through density and vegetation type.
    '''
    
    ### find gamma, hourly melt rate. depends on density 
        #and vegetation type following Kuusisto (1980)

    dd = (9.8e-6 * DENSITY) - 2.39e-3 #daily melt [m water/dayK]
    dd[dd < 1e-4] = 1e-4
    dd[dd > 5.5e-3] = 5.5e-3

    boreal = (5.2e-6 * DENSITY) - 0.7e-3 #calculate daily melt as if all grid squares were Boreal forest [m/dayK]
    boreal[boreal < 1e-4] = 1e-4
    boreal[boreal > 3.5e-3] = 3.5e-3
    
    dd[boreal_mask] = boreal[boreal_mask] #merge dd and boreal

    return dd/24 #melt in [m/hrK]

def reduce_precip(pr, tundraprairie_scaling, tundraprairie_mask, boreal_scaling, boreal_mask):
    '''
    Uniform reduction for two snow class types. Trying
    to capture canopy interception/sublimation effects for
    boreal snowpack. Trying to capture blowing snow sublimation
    loss for tundra and prairie snowpack.
    '''
    pr_copy = pr.copy()
    pr_copy[tundraprairie_mask] = tundraprairie_scaling * pr[tundraprairie_mask]
    pr_copy[boreal_mask] = boreal_scaling * pr[boreal_mask]
   
    return pr_copy

def hourly_melt(weight, mixed_pr, GAMMA, T, PCPN, DEPTH, DENSITY):
    '''
    Compute new snowfall density as a function of air temp
    following Hedstrom and Pomeroy (1998). For temperatures
    greater than 0C, assume linear relationship between snow 
    density and temp following Pomeroy and Gray (1995) Fig. 1
    to max of 200 kg/m^3 (RB 30/09/98)
    
    Args:
        weight: 0.5 used for first and last hour in chunk, 1 used for others
        T: temperature at substep
        PCPN: 1h precipitation (mm water)
        DEPTH: existing snow depth (mm water equivalent)
        DENSITY: existing snow density (kg/m^3)
        GAMMA: melting rate (mm/hr)
    '''
    delt = 3600 #1h [s]
    rhomin = 200 #[kg/m^3]
    rhomax = 550 #[kg/m^3]         
    Tmelt = -1 #[degree C], melt threshold temp
    Tfreeze = 0
    T_switch_upper = 0
    T_switch_lower = 0
    mixed_range = T_switch_upper - T_switch_lower
    
    iopen = np.ones_like(DEPTH)
    icl = np.ones_like(DEPTH)
    
    ### determine precipitation phase at grid squares
    phase = np.where(T <= Tfreeze, 1., 0.) #snow: phase = 1, rain: phase = 0
    if mixed_pr:
        mixed_regime = (T > T_switch_lower) & (T < T_switch_upper)
        phase[mixed_regime] = 1 - (1/mixed_range) * T[mixed_regime]

    SNOW = PCPN * phase #[m water] in one hour
    RAIN = PCPN * (1 - phase) #[m water] in one hour
      
    ### calculate density of new snow based on temperature   
    rhosfall = new_snow_density(T)[SNOW > 0]

    ### change swe units, weight if first or last step
    sfall = 1000 * weight * SNOW[SNOW > 0] #[mm water equivalent]
    
    ### calculate snowpack density through weighted average, with minimum allowed density enforced
    DENSITY[SNOW > 0] = ((rhosfall * sfall) + (DENSITY[SNOW > 0] * DEPTH[SNOW > 0])) / (DEPTH[SNOW > 0] + sfall)
    DENSITY[SNOW > 0] = np.maximum(rhomin, DENSITY[SNOW > 0])

    ### add new snow to depth
    DEPTH[SNOW > 0] += sfall 
    
    ### deal with rain melt
    rain_melt_mask = (RAIN > 0) & (T > Tfreeze) & (DEPTH > 0)
    del_SWE1 = weight * rain_melt(RAIN, T-Tfreeze, rain_melt_mask) #[m water]
    
    ### melt at temperature T    
    temp_melt_mask = (T > Tmelt) & (DEPTH > 0)
    del_SWE2 = weight * temp_melt(T-Tmelt, GAMMA, temp_melt_mask) #[m water]
    
    DEPTH = np.maximum(0., DEPTH - (1000 * del_SWE1) - (1000 * del_SWE2))
    
    ### age snow at T        
    WET_SETTLE = (T >= Tmelt) & (DEPTH > 0)
    COLD_SETTLE = (T < Tmelt) & (DEPTH > 0)
    
    del_DENSITY_warm = warm_snow_aging(DENSITY, DEPTH, delt * weight, WET_SETTLE)
    del_DENSITY_cold = cold_snow_aging(DENSITY, DEPTH, T - Tmelt, icl, COLD_SETTLE)
    
    DENSITY[WET_SETTLE] = np.minimum(rhomax, DENSITY[WET_SETTLE] + del_DENSITY_warm)
    DENSITY[COLD_SETTLE] = np.minimum(rhomax, DENSITY[COLD_SETTLE] + weight * del_DENSITY_cold)
    
    return DEPTH, DENSITY

def Brasnett(mixed_pr, T, pr, OLD, CD, tundraprairie_scaling=0.8, boreal_scaling=0.8):
    '''
    Empirical algorithm to melt snow according to the surface temperature and 
    increase snow depth according to the precipitation that has fallen since 
    the last analysis time.
    
    Args:
        TSFC (floats): ndarray of shape (2, lat, lon) containing the temperature 
            [degree C] at the previous step and current step.
        OLD (float): previous time step snow depth field [m].
        CD (float): previous time step density of the snow pack [kg/m^3]
        PCPN (float): total precipitation [m water] occurring during the time step
    
    #if want to use:
    #Sturm snow classification (icl):
        Water = 0        RHOMIN
        tundra snow = 1    200
        taiga snow = 2     160
        maritime snow = 3  160
        ephemeral snow = 4 180
        prairie snow = 5   140
        alpine snow = 6    120
        ice = 7 (ice caps) 200    
        
    #min density for Sturm snow classes, (starting with water so index matches classification number)
    #rhomin = [1000., 200., 160., 160., 180., 140., 120., 200.] #kg/m^3
    '''       
    rhomin = 200 #[kg/m^3]
    rhomax = 550 #[kg/m^3] 
    
    iopen = np.ones_like(OLD)
    icl = np.ones_like(OLD)
    
    T_switch_upper = 0
    T_switch_lower = 0
    mixed_range = T_switch_upper - T_switch_lower
    
    sdep_max = 6. #[m snow]

    DENSITY = np.maximum(rhomin, np.minimum(rhomax, CD)) #[kg/m^3]
    
    ### no snow and none possible
    no_chance_mask = (T[0,:,:] > T_switch_upper) & (T[-1,:,:] > T_switch_upper) & (OLD <= 0)
    
    boreal_mask = (icl == 2) & (iopen == 0)  
    tundraprairie_mask = (icl == 1) | (icl == 5)
    
    GAMMA = hourly_melt_rate(DENSITY, boreal_mask) #[m/hrK]
    
    pr = reduce_precip(pr, 
                       tundraprairie_scaling, tundraprairie_mask,
                       boreal_scaling, boreal_mask)
        
    NEW = OLD * DENSITY #[mm water]
    
    ### beyond this point, NEW and DENSITY will be updated for each hour 
        #based on the temperature and precipitation
    NEW, DENSITY = hourly_melt(0.5, mixed_pr, GAMMA, T[0,:,:], pr, NEW, DENSITY)
    for i in range(1, np.shape(T)[0]-1):
        NEW, DENSITY = hourly_melt(1, mixed_pr, GAMMA, T[i,:,:], pr, NEW, DENSITY)
    NEW, DENSITY = hourly_melt(0.5, mixed_pr, GAMMA, T[-1,:,:], pr, NEW, DENSITY)
    
    ### output after chunk
    NEW = (NEW / DENSITY) #[m snow], total depth of snow layer
    NEW = np.minimum(NEW, sdep_max) #depth does not exceed 6m
    
    NEW[no_chance_mask] = 0.     
    DENSITY[no_chance_mask] = rhomin #[kg/m^3], snowpack density
    
    swe = NEW * DENSITY #[mm water], snowpack water equivalent
    
    return NEW, DENSITY, swe
