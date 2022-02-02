**Contact: Aleksandra Elias Chereque // aleksandra.eliaschereque@mail.utoronto.ca**

## Brown Temperature Index Model (B-TIM)

This simple model uses (sub-daily or daily) temperature and precipitation inputs and generates a daily record of snow depth and snow density on the same grid as the input data. Data start on August 1 and end on July 31, and each hydrological year initializes from 0mm SWE everywhere.


The model captures a single snow layer which evolves under the following processes:
* New snowfall -> Hedstrom/Pomeroy (1998) relationship for new snow density; snow layer density updated as weighted average
* Rain-induced melt 
* Warm-wet snow settling/densification
* Cold snow settling/densification -> follows Anderson (1978)
* Melt due to air-temperatures above -1C -> melt rate parameterization follows Kuusisto (1980)

## Running the model

The python environment can be installed from the environment.yml file with conda:

```
conda env create -f environment.yml
```

BTIM.py contains the main model script, which can be run in the environment for the snow year (Aug YYYY, July YYYY+1) as follows:

```
conda activate env
python BTIM.py YYYY
```

### Working example:
* Run one month with ERA5 forcing that is provided by executing:
```
python BTIM.py 2020
```
* The model run will end with a FileNotFound error for this sample - only August forcing files are found in 'forcing/'
* One output file should have been generated, to be found in 'output/'

### Running with new forcing

To adjust file I/O details, you must make changes in the config.py file and create a forcing-specific file (such as ERA5_BTIM.py).

1. The helper file should contain all of the variables described below as well as four functions.
   **Variables**:
   
   |Name|Description|Example|Options(where relevant)|
   |---|---|---|---|
   |forcing|name of the forcing dataset, typically found in the forcing filenames|'ERA5'|--|
   |data_loc|direct path to forcing data|'/home/forcing_data/'|--|
   |pr_freq|number of time steps per day with new precipitation data|24|1, 2, 3, 4, 6, 8, 12, 24|
   |t2m_freq|number of time steps per day with new temperature data|24|1, 2, 3, 4, 6, 8, 12, 24|
   |pr_in_mmhr|boolean confirming the units of input preciptiation|True|True or False|
   |lat_increasing|for gridded data, confirms if the first row to be read is northmost or southmost latitude|True|True or False|
   |lon_0_360|confirms if longitudes range from -180 to 180 or from 0 to 360|True|True or False|
   |tp_name|name of precipitation variable in input file|'tp'|--|
   |t2m_name|name of temperature variable in input file|'t2m'|--|
   |latitude_name|name of latitude variable in input file|'lat'|--|
   |longitude_name|name of longitude variable in input file|'lon'|--|
   
   **Functions**: 
   
      1. **prepare_filenames** must take in a snow season, e.g. (YYYY,YYYY+1), and return two lists of filenames spanning from August 1, YYYY to July 31, YYYY+1. The first list contains the temperature forcing filenames, and the second list contains the preciptiation forcing filenames.
      2. **read_month** must take in a filename (or multiple filenames) for each variable along with a month (1=Jan). It should load the relevant data from the files. I have loaded into netCDF4 Datasets or into xarray Datasets. The function should return the latitude and longitude grid associated with the data and two datasets containing the forcing data (one for temperature, one for precipitation).
      3. **read_day** must take in datasets as returned by read_month and access a particular time step, indexing 0 as the first time step of the first day of a specified month.
      4. **standardize_precip** should be able to take in (precipitation) data and rescale so that the units are in metres per temperature step.
   
2. With all this in place, editing the config.py file is straightforward for a given experiment. Load in the helper file you just made in the first line (e.g. 'from ERA5_BTIM import \*') and make changes to any varaibles, defined below.

   **Variables:**
   
   |Name|Description|Example|Options|
   |---|---|---|---|
   |unique_ID|a unique name for this experiment/run of the model|'ERA5_rescaled_precip'|--|
   |mixed_pr|turns on and off mixed precipitation for temperatures between 0-2C|True|True or False|
   |latminmax|specifies the latitude range to model|\[40,90\]|Tuple with values between -90 and 90|
   |lonminmax|specifies the longitude range to model|\[0,360\]|Tuple with valid longitude values, wraparound supported (e.g. \[270, 25\]|
   |leapdays|confirms whether the forcing data calendar has leap days|True|True or False|
   |output_loc|path to directory where files should be saved|'./'|--|


## References
Original mention of the algorithm:
* Brasnett, B. (1999). A Global Analysis of Snow Depth for Numerical Weather Prediction, Journal of Applied Meteorology, 38(6), 726-740. Retrieved Jan 31, 2022, from https://journals.ametsoc.org/view/journals/apme/38/6/1520-0450_1999_038_0726_agaosd_2.0.co_2.xml

Update:
* Ross D. Brown, Bruce Brasnett & David Robinson (2003) Gridded North American monthly snow depth and snow water equivalent for GCM evaluation, Atmosphere-Ocean, 41:1, 1-14, DOI: 10.3137/ao.410101

BTIM + snow depth observation OI used:
* Brown, R. D., & Mote, P. W. (2009). The Response of Northern Hemisphere Snow Cover to a Changing Climate, Journal of Climate, 22(8), 2124-2145. Retrieved Jan 31, 2022, from https://journals.ametsoc.org/view/journals/clim/22/8/2008jcli2665.1.xml

The most recent user guide for the CMC snow depth dataset, which uses this model for the first guess field:
* https://nsidc.org/sites/nsidc.org/files/NSIDC-0447-V001-UserGuide_0.pdf
