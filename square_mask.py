import numpy as np

def square_mask(full_lat, full_lon, latminmax, lonminmax, lon_0_360):
    
    lat_mask = (full_lat >= latminmax[0]) & (full_lat <= latminmax[1])
    
    if lon_0_360:
        lon_mask = (full_lon >= np.min(lonminmax)) & (full_lon <= np.max(lonminmax))
    
    else:
        if lonminmax[0] < lonminmax[1]:
            if lonminmax[0] < 180:
                if lonminmax[1] >= 180:
                    lon_mask = (full_lon >= lonminmax[0]) | (full_lon <= (lonminmax[1]+180)%360-180)
                else:
                    lon_mask = (full_lon >= lonminmax[0]) | (full_lon <= lonminmax[1])
            else:
                lon_mask = (full_lon >= (lonminmax[0]+180)%360-180) | (full_lon <= (lonminmax[1]+180)%360-180)                   
        else:
            if lonminmax[1] < 180:
                lon_mask = (full_lon >= lonminmax[0]) | (full_lon <= lonminmax[1])
            else:
                lon_mask = (full_lon >= lonminmax[0]) | (full_lon <= 180) | (full_lon <= (lonminmax[1]+180)%360-180)

    return lat_mask, lon_mask
