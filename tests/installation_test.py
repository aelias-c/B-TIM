import xarray as xr

data = xr.open_dataset('test_result.nc')

for v in data.data_vars:
    if data[v].mean() < 1e-6:
        print("TEST PASSED! var: "+v)
    else:
        print("TEST FAILED! var: "+v)
    