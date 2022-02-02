import numpy as np

def snow_aging(DENSITY, DEPTH, T, icl, Tmelt, delt, weight, rhomax): 
    TDD = T - Tmelt
    # two regimes: warm/wet settling, cold settling
    WET_SETTLE = (T >= Tmelt) & (DEPTH > 0)
    COLD_SETTLE = (T < Tmelt) & (DEPTH > 0.)                           

    #warm, wet snow
    sdepcm = 100 * (DEPTH[WET_SETTLE] / DENSITY[WET_SETTLE])  #[cm]
    denmax = 700. - ((20470. / sdepcm) * (1 - np.exp(-sdepcm / 67.3)))
    den_diff = denmax - DENSITY[WET_SETTLE]

    TIMFAC = np.exp(np.log(den_diff[den_diff > 0.1] / 200.) - (2.778e-6 * delt * weight))

    DENSITY[WET_SETTLE][den_diff> 0.1] = denmax[den_diff > 0.1] - 200. * TIMFAC
    DENSITY[WET_SETTLE] = np.minimum(rhomax, DENSITY[WET_SETTLE])

    #cold snow aging
    C2 = np.where((icl == 2) | (icl == 6), -28, -21)[COLD_SETTLE]
    den_gcm = DENSITY[COLD_SETTLE] / 1000
    del_den = 0.02 * np.exp(0.08 * TDD[COLD_SETTLE]) * 0.6 * DEPTH[COLD_SETTLE] * np.exp(C2 * den_gcm) / 10
    dgain = del_den * 1000

    DENSITY[COLD_SETTLE] = DENSITY[COLD_SETTLE] + (dgain * weight)
    DENSITY[COLD_SETTLE] = np.minimum(rhomax, DENSITY[COLD_SETTLE])

    return DENSITY                          

def temp_melt(DEPTH, T, GAMMA, Tmelt, weight):
    TDD = T - Tmelt
    TEMP_MELT = (TDD > 0) & (DEPTH > 0)

    DEPTH[TEMP_MELT] = DEPTH[TEMP_MELT] - (weight * GAMMA[TEMP_MELT] * TDD[TEMP_MELT])
    
    return DEPTH

def rain_melt(RAIN, weight, T, DEPTH, rain_melt_mask, delt):
    rhow = 1000 #[kg/m^3], density of water
    Cw = 4.18e3 #[J], specific heat of water
    rhoice = 917 #[kg/m^3], density of ice
    Lf = 0.334e6 #[J/kg], latent heat of fusion of water

    prain = np.where(rain_melt_mask, RAIN * 1000 / delt, 0.)
    rmelt = (rhow * Cw * prain * T) / (Lf * rhoice) #[mm/s]

    DEPTH = DEPTH - (weight * delt * rmelt)
    DEPTH = np.maximum(0., DEPTH)

    return DEPTH


def new_snow_density(T):
    '''
    New snow density follows Hedstrom and Pomeroy (1998).
    '''

    rhosfall_cold = 67.9 + 51.3 * np.exp(T / 2.6)
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

    dd = ((0.0196 / 2) * DENSITY) - 2.39 #daily melt [mm/day]
    dd[dd < 0.1] = 0.1
    dd[dd > 5.5] = 5.5

    boreal = ((0.0104 / 2) * DENSITY) - 0.70 #calculate daily melt as if all grid squares were Boreal forest [mm/day]
    boreal[boreal < 0.1] = 0.1 
    boreal[boreal > 3.5] = 3.5


    dd[boreal_mask] = boreal[boreal_mask] #merge dd and boreal

    return dd/24 #melt in [mm/hr]

def reduce_precip(pr, boreal_scaling, boreal_mask, tundraprairie_scaling, tundraprairie_mask):
    '''
    Uniform reduction for two snow class types. Trying
    to capture canopy interception/sublimation effects for
    boreal snowpack. Trying to capture blowing snow sublimation
    loss for tundra and prairie snowpack.
    '''
    pr[boreal_mask] = boreal_scaling * pr[boreal_mask]
    pr[tundraprairie_mask] = tundraprairie_mask * pr[tundraprairie_mask]

    return pr

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
    
    iopen = 1
    icl = np.ones_like(DEPTH)
    
    ### determine precipitation phase at grid squares
    phase = np.where(T <= 0, 1., 0.) #snow: phase = 1, rain: phase = 0
    if mixed_pr:
        mixed_regime = (T > 0) & (T < 2)
        phase[mixed_regime] = 1 - 0.5 * T[mixed_regime]

    SNOW = PCPN * phase #[m water] in one hour
    RAIN = PCPN * (1 - phase) #[m water] in one hour
      
    ### calculate density of new snow based on temperature   
    rhosfall = new_snow_density(T)[SNOW > 0]

    ### change swe units, weight if first or last step
    sfall = weight * 1000 * SNOW[SNOW > 0] #[mm water equivalent]

    ### add new snow to depth
    DEPTH[SNOW > 0] += sfall 
    
    ### calculate snowpack density through weighted average, with minimum allowed density enforced
    DENSITY[SNOW > 0] = ((rhosfall * sfall) + DENSITY[SNOW > 0] * (DEPTH[SNOW > 0] - sfall)) / DEPTH[SNOW > 0] 
    DENSITY[SNOW > 0] = np.maximum(rhomin, DENSITY[SNOW > 0])

    ### deal with rain melt
    rain_melt_mask = (RAIN > 0) & (T > 0) & (DEPTH > 0)
    DEPTH = rain_melt(RAIN, weight, T, DEPTH, rain_melt_mask, delt)
                 
    ### melt at temperature T       
    DEPTH = np.maximum(0., temp_melt(DEPTH, T, GAMMA, Tmelt, weight))
    
    ### age snow at T                        
    DENSITY = snow_aging(DENSITY, DEPTH, T, icl, Tmelt, delt, weight, rhomax)
    
    return DEPTH, DENSITY

def Brasnett(hours_per_chunk, mixed_pr, T, pr, OLD, CD, boreal_scaling=0.8, tundraprairie_scaling=0.8):
    '''
    Empirical algorithm to melt snow according to the surface temperature and 
    increase snow depth according to the precipitation that has fallen since 
    the last analysis time.
    
    Args:
        TSFC (floats): ndarray of shape (2, lat, lon) containing the temperature 
            [degree C] at the previous step and current step.
        OLD (float): previous time step snow depth field [cm].
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
    iopen = 1
    icl = 1
    sdep_max = 600 #[cm snow]

    DENSITY = np.maximum(rhomin, np.minimum(rhomax, CD)) #[kg/m^3]
    
    ### no snow and none possible
    no_chance_mask = (T[0,:,:] > 2) & (T[-1,:,:] > 2) & (OLD <= 0)
    
    boreal_mask = (icl == 2) & (iopen == 0)  
    tundraprairie_mask = (icl == 1) | (icl == 5)
    
    GAMMA = hourly_melt_rate(DENSITY, boreal_mask)
    
    pr = reduce_precip(pr, boreal_scaling, boreal_mask, 
                      tundraprairie_scaling, tundraprairie_mask)
        
    NEW = OLD * DENSITY * 0.01 #[mm water]
    
    ### beyond this point, NEW and DENSITY will be updated for each hour 
        #based on the temperature and precipitation
    NEW, DENSITY = hourly_melt(0.5, mixed_pr, GAMMA, T[0,:,:], pr, NEW, DENSITY)
    for i in range(1, hours_per_chunk):
        NEW, DENSITY = hourly_melt(1, mixed_pr, GAMMA, T[i,:,:], pr, NEW, DENSITY)
    NEW, DENSITY = hourly_melt(0.5, mixed_pr, GAMMA, T[-1,:,:], pr, NEW, DENSITY)
    
    ### output after chunk
    NEW = 100 * (NEW / DENSITY) #[cm snow], total depth of snow layer
    NEW = np.minimum(NEW, sdep_max) #depth does not exceed 6m
    
    NEW[no_chance_mask] = 0.     
    DENSITY[no_chance_mask] = rhomin #[kg/m^3], snowpack density
    
    swe = NEW * 0.01 * DENSITY #[mm water], snowpack water equivalent
    
    return NEW, DENSITY, swe
