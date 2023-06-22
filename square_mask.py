def square_mask(full_lat, full_lon, latminmax, lonminmax):
    
    # mask for latitudes
    lat_mask = (full_lat >= latminmax[0]) & (full_lat <= latminmax[1])
    
    # mask for longitudes
    lon_mask = (full_lon >= 0) | (full_lon <= 0) 

    return lat_mask, lon_mask
                
       
