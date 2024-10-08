import numpy as np

constants = {
    'rhow': 1000., #[kg/m3], density of water
    'rhoice': 917., #[kg/m3], density of ice
    'Cw' : 4.18e3, #[J], specific heat of water
    'Lf' : 0.334e6, #[J/kg], latent heat of fusion of water\
}

params = {   
    'rhomin': 200., #[kg/m^3], min allowed density
    'rhomax': 550., #[kg/m^3], max allowed density     
    'sdep_max' : 6., #[m snow], max allowed depth
    'Tmelt': -1., #[degree C], melt threshold temp
    'Tfreeze': 0. #[degree C]
}

def warm_snow_aging(DENSITY, DEPTH, timestep, warm_settle_mask):
    '''
    Density change over timestep when temperature exceeds T_melt.
    
    Args:
        DENSITY (ndarray): current snow density [kg/m3]
        DEPTH (ndarray): current snow depth [m snow]
        timestep (float): duration of warm snow aging [s]
        warm_settle_mask (ndarray): boolean mask of the same size 
            as DENSITY and SWE denoting where warm snow aging should 
            occur.

    Returns:
        del_den_warm (ndarray): density increment valid wherever 
            wet_settle_mask has value 1.
    '''

    Wmax, W1, W2 = 700., 204.70, 0.673
    a = 2.778e-6
    
    denmax = Wmax - ((W1 / DEPTH[warm_settle_mask]) * (1 - np.exp(-DEPTH[warm_settle_mask] / W2)))
    den_diff = np.atleast_1d(denmax - DENSITY[warm_settle_mask])

    del_den_warm = np.zeros_like(denmax)
    del_den_warm[den_diff > 0.1]  =  den_diff[den_diff > 0.1] * (1 - np.exp(-a * timestep))
    del_den_warm = np.atleast_1d(del_den_warm)

    return del_den_warm

def cold_snow_aging(DENSITY, DEPTH, TDD, icl, cold_settle_mask):
    '''
    Density change over timestep when temperature is below T_melt.
    
    Args:
        DENSITY (ndarray): current snow density [kg/m3]
        DEPTH (ndarray): current snow depth (m snow)
        TDD (ndarray): current temperature minus T_melt
        icl (ndarray): land cover classification of same 
            size as above arrays
        cold_settle_mask (ndarray): boolean mask of the same size 
            as DENSITY and SWE denoting where cold snow aging should 
            occur.

    Returns:
        del_den_cold (ndarray): density increment valid wherever 
            cold_settle_mask has value 1.
    '''

    C1 = 2.
    C2 = np.where((icl == 2) | (icl == 6), 28./1000., 21./1000.)
    C3 = 0.08
    B1 = 0.6
    
    del_den_cold = C1 * (B1 * DENSITY[cold_settle_mask] * DEPTH[cold_settle_mask]) * np.exp(C3 * TDD[cold_settle_mask]) * np.exp(-C2[cold_settle_mask] * DENSITY[cold_settle_mask]) #[kg/m3]
    
    return del_den_cold                          

def temp_melt(TDD, DENSITY, hG, temp_melt_mask):
    '''
    SWE change due to melting (occurs when temperature exceeds T_melt).
    
    Args:
        TDD (ndarray): current temperature minus T_melt
        DENSITY (ndarray): current snow density [kg/m3]
        hG (ndarray): melt rate (mm w.e./hrK), array of same size as TDD
        temp_melt_mask (ndarray): boolean mask of the same size as TDD, hG
            denoting where temperature melt should occur.

    Returns:
        del_DEPTH (ndarray): array containing snow depth increments (m snow)
    '''

    del_DEPTH = np.where(temp_melt_mask, TDD*hG/DENSITY, 0.)
    
    return del_DEPTH

def rain_melt(RAIN, DENSITY, T_diff, rain_melt_mask):    
    '''
    Snow depth change due to melting (occurs when rain falls on snow).
    
    Args:
        RAIN (ndarray): liquid precipitation falling during time step (m w.e.)
        DENSITY (ndarray): current snow density [kg/m3]
        T_diff (ndarray): current temperature minus T_freeze (=0C).
        rain_melt_mask (ndarray): boolean mask of the same size as RAIN, T_diff
            denoting where rain melt should occur.

    Returns:
        del_DEPTH (ndarray): array containing snow depth increments (m snow)
    '''
    
    ### RAIN has units [m water] equiv. [1 m3 water per m2]
    mass_water_per_m2 = np.where(rain_melt_mask, constants['rhow'] * RAIN, 0.) #[kg/m2] 
    heat_from_rain = np.where(rain_melt_mask, mass_water_per_m2 * constants['Cw'] * T_diff, 0.) #[J/m2]

    del_DEPTH = np.where(rain_melt_mask, heat_from_rain / (constants['Lf'] * DENSITY), 0.) #[m snow]

    return del_DEPTH


def new_snow_density(hT):
    '''
    Density of new snow, following Hedstrom and Pomeroy (1998).
    
    Args:
        hT (ndarray): current temperature in degreesC

    Returns: 
        rhosfall (ndarray): density [kg/m3] of fresh snow, array of 
            same size as hT.
    '''
    
    rhosfall_cold = 67.9 + 51.3 * np.exp(hT / 2.6)
    
    ### rhosfall_warm will only be relevant ###
    ###     if mixed precip is active       ###
    rhosfall_warm =  np.minimum(119.2 + (20 * hT), 200.)
  
    rhosfall = np.where(hT <= 0, rhosfall_cold, rhosfall_warm)
    ### ----------------------------------- ###
    
    return np.atleast_1d(rhosfall)

def daily_melt_rate(DENSITY, boreal_mask):
    ''' 
    Daily melt rate for snow based on existing snow density 
    and vegetation type. Based on Kuusito (1980).

    Args:
        DENSITY (ndarray): current density of snow [kg/m3].
        boreal_mask (ndarray): boolean array of the same size as DENSITY
            denoting regions which are classified as boreal forest land
            type.

    Returns:
        dd (ndarray): hourly melt rate [mm w.e. / hr.day] of snowpack.
    '''
    
    ### find gamma, hourly melt rate. depends on density 
        #and vegetation type following Kuusisto (1980)

    dd = (9.8e-3 * DENSITY) - 2.39 #daily melt [mm w.e./dayK]
    dd[dd < 0.1] = 0.1
    dd[dd > 5.5] = 5.5

    ### AEC 2023: we're treating all points as Tundra  ###
    ### so this section is basically ignored           ###
    boreal = (5.2e-3 * DENSITY) - 0.7 #[mm w.e./dayK]
    boreal[boreal < 0.1] = 0.1
    boreal[boreal > 3.5] = 3.5
    ### ---------------------------------------------- ###
    
    dd = np.where(boreal_mask, boreal, dd) #merge dd and boreal

    return dd #melt in [mm w.e./ hr day]

def reduce_precip(hP, tundraprairie_scaling, tundraprairie_mask, boreal_scaling, boreal_mask):
    '''
    Uniform reduction of precipitation for two snow class types. 
    Implemented to capture canopy interception/sublimation effects for
    boreal snowpack and blowing snow sublimation loss for tundra and 
    prairie snowpack.

    Args:
        hP (ndarray): precipitation for current hour (mm w.e.)
        tundraprairie_scaling (float): value between 0 and 1 by which to 
            multiplicatively scale precipitation over tundra/prairie snowpack
        tundraprairie_mask (ndarray): boolean mask of the same size as hP 
            denoting regions classified as tundra or prairie
        boreal_scaling (float): value between 0 and 1 by which to 
            multiplicatively scale precipitation over boreal snowpack
        boreal_mask (ndarray): boolean mask of the same size as hP 
            denoting regions classified as boreal forest

    Returns:
        hP (ndarray): scaled precipitation
    '''
    
    hP = np.where(tundraprairie_mask, hP * tundraprairie_scaling, hP)
    hP = np.where(boreal_mask, hP * boreal_scaling, hP)
   
    return hP

def hour_step(weight, mixed_pr_range, hG, hT, hP, DENSITY, DEPTH, debug_mode=False):
    '''
    Update snow density given preceipitation and temperature
    during time step.
    
    If mixed precipitation is allowed, assume linear relationship between snow 
    density and temp above 0C following Pomeroy and Gray (1995) Fig. 1
    (RB 30/09/98)
    
    Args:
        weight (float): value between 0 and 1. Used to center
            calculations around the top of the hour. 0.5 used for first
            and last forcing values in one forcing time step. 1 used 
            for full hours within a forcing time step.
        mixed_pr_range (tuple): lower and upper threshold (degreeC) temperatures
            for mixed precipitation.
        hG: hourly melting rate (mm w.e.)
        hT (ndarray): hourly mean temperature
        hP (ndarray): hourly precipitation (mm w.e.)
        DENSITY: existing snow density (kg/m^3)
        DEPTH: existing snow depth (m snow)

    Returns:
        DEPTH (ndarray): updated depth, array of same size as hT and hP
        DENSITY (ndarray): updated density, array of same size as hT and hP
    '''

    if debug_mode:
        print('hour step init (T, P, rho, D):', hT, hP, DENSITY, DEPTH)
    
    T_switch_upper, T_switch_lower = mixed_pr_range
    mixed_range = T_switch_upper - T_switch_lower
    
    iopen = np.ones_like(DENSITY)
    icl = np.ones_like(DENSITY)
    
    ### determine precipitation phase at grid squares
    phase = np.where(hT <= params['Tfreeze'], 1., 0.) #snow: phase = 1, rain: phase = 0
    if T_switch_lower != T_switch_upper:
        mixed_regime = (hT > T_switch_lower) & (hT < T_switch_upper)
        phase[mixed_regime] = 1 - (1/mixed_range) * hT[mixed_regime]

    SNOW = np.atleast_1d(hP * phase) #[m water] in one hour
    RAIN = np.atleast_1d(hP * (1 - phase)) #[m water] in one hour

    ### calculate density of new snow based on temperature   
    rhosfall = np.where(SNOW > 0, new_snow_density(hT), 0.)

    ### change swe units, weight if first or last step
    swefall = np.atleast_1d(weight * constants['rhow'] *  np.where(SNOW > 0, SNOW, 0.)) #[mm water equivalent]

    SWE0 = np.atleast_1d(DEPTH * DENSITY)
    
    ### calculate snowpack density through weighted average, with minimum allowed density enforced
    DENSITY = np.where((SNOW > 0), 
                       np.divide(np.atleast_1d(rhosfall * swefall) + (DENSITY * SWE0),
                                 (SWE0 + swefall),
                                 out=np.zeros_like(DENSITY),
                                 where=(SWE0+swefall!=0)),
                       DENSITY)
    DENSITY = np.minimum(params['rhomax'], np.maximum(params['rhomin'], DENSITY))

    ### add new snow to depth
    SWE0 = SWE0 + swefall
    DEPTH = SWE0 / DENSITY

    ### shape
    DEPTH = np.atleast_1d(DEPTH)

    if debug_mode:
        print('New rho and D:', DENSITY, DEPTH)

    ### deal with rain melt
    rain_melt_mask = np.atleast_1d( (RAIN > 0) & (DEPTH > 0) & (hT > params['Tfreeze']) )

    del_DEPTH1 = -rain_melt(np.atleast_1d(RAIN), np.atleast_1d(DENSITY), np.atleast_1d(hT-params['Tfreeze']), rain_melt_mask) #[m snow]

    if debug_mode:
        print('Delta rain melt:', del_DEPTH1)

    ### melt at temperature T    
    temp_melt_mask = np.atleast_1d((hT > params['Tmelt']) & (DEPTH > 0))
    del_DEPTH2 = -temp_melt(np.atleast_1d(hT-params['Tmelt']), DENSITY, np.atleast_1d(hG), temp_melt_mask) #[m snow]

    if debug_mode:
        print('DENSITY, GAMMA:', DENSITY, hG)
        print('Delta melt:', del_DEPTH2)
    
    DEPTH = np.atleast_1d(np.maximum(0., DEPTH + weight * del_DEPTH1 + weight * del_DEPTH2))

    if debug_mode:
        print('Dstarstar:', DEPTH)
    
    ### age snow at T   
    SWEf = np.atleast_1d(DEPTH * DENSITY)

    warm_settle_mask = np.atleast_1d( (hT >= params['Tmelt']) & (DEPTH > 0.) )
    cold_settle_mask = np.atleast_1d( (hT < params['Tmelt']) & (DEPTH > 0.) )

    btim_tdelt = 3600 #1h [s]
    del_DENSITY_warm = warm_snow_aging(DENSITY, DEPTH, btim_tdelt * weight, warm_settle_mask)
    
    del_DENSITY_cold = weight * cold_snow_aging(DENSITY, DEPTH, np.atleast_1d(hT - params['Tmelt']), icl, cold_settle_mask)

    del_DENSITY = np.zeros_like(DENSITY)

    del_DENSITY[warm_settle_mask] = del_DENSITY_warm    
    del_DENSITY[cold_settle_mask] = del_DENSITY_cold    
    DENSITY = np.maximum(params['rhomin'], np.minimum(params['rhomax'], DENSITY + del_DENSITY))

    if debug_mode:
        print('Del rho, densification', del_DENSITY)

    DEPTH = np.atleast_1d(SWEf / DENSITY) # conserve water
   
    if debug_mode:
        print('Hour end rho, D:', DENSITY, DEPTH)
    
    return DENSITY, DEPTH

def Brasnett(mixed_pr_range, HOURLY_T, HOURLY_PRECIP, SNOW_DEPTH, SNOW_DENSITY, tundraprairie_scaling=0.8, boreal_scaling=0.8, debug_mode=False):
    '''
    Empirical algorithm to melt snow according to the surface temperature and 
    increase snow depth according to the precipitation that has fallen since 
    the last analysis time.
    
    Args:
        HOURLY_T (floats): ndarray of shape (nhours+1, lat, lon) containing the temperature 
            [degree C] every hour during one forcing data time step, which may be longer
            than one hour.
        HOURLY_PRECIP (float): total precipitation [m water] occurring per hour during time step
        SNOW_DEPTH (float): snow depth field [m] at the beginning of the time step.
        SNOW_DENSITY (float): density field at the beginning of the time step [kg/m^3]
    '''     
    
    iopen = np.ones_like(SNOW_DEPTH)
    icl = np.ones_like(SNOW_DEPTH)
    
    T_switch_upper, T_switch_lower = mixed_pr_range
    mixed_range = T_switch_upper - T_switch_lower
    
    SNOW_DENSITY = np.maximum(params['rhomin'], np.minimum(params['rhomax'], SNOW_DENSITY)) #[kg/m^3]
    
    ### no snow and none possible
    no_chance_mask = (HOURLY_T[0,...] > T_switch_upper) & (HOURLY_T[-1,...] > T_switch_upper) & (SNOW_DEPTH <= 0)

    ### calculate melt rates
    boreal_mask = (icl == 2) & (iopen == 0)  
    tundraprairie_mask = (icl == 1) | (icl == 5)

    HOURLY_GAMMA = daily_melt_rate(SNOW_DENSITY, boreal_mask) / 24 #[mm/hrK]

    ### adjust precipitation
    HOURLY_PRECIP = reduce_precip(HOURLY_PRECIP, 
                                  tundraprairie_scaling, tundraprairie_mask,
                                  boreal_scaling, boreal_mask)
    
    ### beyond this point, SNOW_DEPTH and SNOW_DENSITY will be updated for each hour in the 
        #model time step based on temperature and precipitation
    SNOW_DENSITY, SNOW_DEPTH = hour_step(0.5, mixed_pr_range, HOURLY_GAMMA, 
                                         np.atleast_1d(HOURLY_T[0,...]), HOURLY_PRECIP, SNOW_DENSITY, SNOW_DEPTH,
                                         debug_mode=debug_mode)
    for i in range(1, np.shape(HOURLY_T)[0]-1):
        SNOW_DENSITY, SNOW_DEPTH = hour_step(1, mixed_pr_range, HOURLY_GAMMA, 
                                             np.atleast_1d(HOURLY_T[i,...]), HOURLY_PRECIP, SNOW_DENSITY, SNOW_DEPTH,
                                             debug_mode=debug_mode)
    SNOW_DENSITY, SNOW_DEPTH = hour_step(0.5, mixed_pr_range, HOURLY_GAMMA, 
                                         np.atleast_1d(HOURLY_T[-1,...]), HOURLY_PRECIP, SNOW_DENSITY, SNOW_DEPTH, 
                                         debug_mode=debug_mode)
    
    ### save final value after model time step
    SNOW_DEPTH = np.minimum(SNOW_DEPTH, params['sdep_max']) #depth does not exceed 6m
    
    SNOW_DEPTH = np.where(no_chance_mask, 0., SNOW_DEPTH)     
    SNOW_DENSITY = np.where(no_chance_mask, params['rhomin'], SNOW_DENSITY) #[kg/m^3], snowpack density
    
    NEW_SWE = SNOW_DEPTH * SNOW_DENSITY #[mm water], snowpack water equivalent
    
    return SNOW_DEPTH, SNOW_DENSITY, NEW_SWE
