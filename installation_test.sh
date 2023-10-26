#!/bin/bash/

python BTIM.py 2019 ERA5
ncdiff -O output/ERA5_forced_swe_Aug_2019_2020.nc  output/out_test.nc test_result.nc
python tests/installation_test.py
