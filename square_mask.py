def square_mask(full_lat, full_lon, latminmax):
    
    # mask for latitudes
    lat_mask = (full_lat >= latminmax[0]) & (full_lat <= latminmax[1])
    
    return lat_mask
                
       
