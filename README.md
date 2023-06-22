**Contact: Aleksandra Elias Chereque // aleksandra.eliaschereque@mail.utoronto.ca**

## Brown Temperature Index Model (B-TIM)

This simple distributed model uses (sub-daily or daily frequency) temperature and precipitation inputs and generates a daily record of snow depth and snow density on the same grid as the input data. Each water year starts on August 1 and ends on July 31, with 0m snow depth everywhere on August 1. The current configuration accepts square grids or single points only.

The model captures a single snow layer which evolves under the following processes:
* New snowfall -> Hedstrom/Pomeroy (1998) relationship for new snow density; snow layer density updated as weighted average
* Rain-induced melt 
* Warm-wet snow settling/densification OR Cold snow settling/densification -> follows Anderson (1978)
* Melt due to air-temperatures above -1C -> melt rate parameterization follows Kuusisto (1980)

## Running the model

The python environment can be installed from the environment.yml file with conda:

```
conda env create -f environment.yml
```

BTIM.py contains the main model script, which can be run in the environment for the snow year (Aug YYYY, July YYYY+1) and forcing F as follows:

```
conda activate env
python BTIM.py YYYY F
```

### Working example:
* Run one month with the provided ERA5 forcing by executing:

```
python BTIM.py 2020
```

* _The model run will end with a FileNotFound error for this sample - this is because there are only August forcing files in 'forcing/'_
* One output file should have been generated, which will be found in 'output/'

### Running with new forcing

For other applications, it will be necessary to format the input files and configuration details as follows:

1. The helper file should contain all of the variables described below as well as four functions.
   **Variables in CONFIG.py**:
   
   |Name|Description|Example|Options(where relevant)|
   |---|---|---|---|
   |data_loc|direct path to forcing data|'/home/forcing_data/'|--|
   |mixed_pr|boolean flag to include some frozen precipitation between 0-2C|True|True or False|
   |latminmax|lower and upper latitude bounds (inclusive) for a rectangular subset of the globe, tuple of size (2,)|\[10, 90\]|--|
   |lonminmax|lower and upper longitude bounds (from 0-360, inclusive) for a rectangular subset of the globe, tuple of size (2,)|\[0, 360\]|--|
   |leapdays|boolean flag, does the forcing data have leap days?|True|True or False|
   |output_loc|path to directory where files should be saved|'./'|--|
    |unique_ID|a unique name for this experiment/run of the model|'ERA5_rescaled_precip'|--|
   
    **Settings in utils.py**:
   
   * Add forcing name and data to each dictionary following the pattern provided (t2m_freq, tp_freq, latname, lonname, precipname, tempname)
   * tp_freq: number of time steps per day with new precipitation data (possible values: 1, 2, 3, 4, 6, 8, 12, 24)
   * t2m_freq: number of time steps per day with new temperature data (possible values: 1, 2, 3, 4, 6, 8, 12, 24)
   * precipname: name of precipitation variable in input file
   * tempname: name of temperature variable in input file
   * latname: name of latitude variable in input file
   * lonname: name of longitude variable in input file
 
   **Functions in utils.py**: 
   
      1. **prepare_filenames** if your forcing files follow a different naming convention, update this function. It must take in a snow season, e.g. (YYYY,YYYY+1), and return two lists of filenames spanning from August 1, YYYY to July 31, YYYY+1. The first list contains the temperature forcing filenames, and the second list contains the preciptiation forcing filenames.
      2. **read_month** must take in a filename (or multiple filenames) for each forcing variable along with a month (1=Jan) and year. It should load the relevant data from the files. I have loaded into netCDF4 Datasets or into xarray Datasets. The function should also return the latitude and longitude grid associated with the data and two datasets containing the forcing data (one for temperature, one for precipitation).
      3. **read_day** must take in datasets as returned by read_month and access a particular time step, indexing 0 as the first time step of the first day of a specified month. The data should be returned as an array.
      4. **standardize_precip** should be able to take in (precipitation) data and rescale so that the units are in metres per temperature step.
      5. **monthly_out_name** if desired, adjust this to change how the files are named.

## References
Original mention of the algorithm:
* Brasnett, B. (1999). A Global Analysis of Snow Depth for Numerical Weather Prediction, Journal of Applied Meteorology, 38(6), 726-740. Retrieved Jan 31, 2022, from https://journals.ametsoc.org/view/journals/apme/38/6/1520-0450_1999_038_0726_agaosd_2.0.co_2.xml

Update:
* Ross D. Brown, Bruce Brasnett & David Robinson (2003) Gridded North American monthly snow depth and snow water equivalent for GCM evaluation, Atmosphere-Ocean, 41:1, 1-14, DOI: 10.3137/ao.410101

BTIM + snow depth observation OI used:
* Brown, R. D., & Mote, P. W. (2009). The Response of Northern Hemisphere Snow Cover to a Changing Climate, Journal of Climate, 22(8), 2124-2145. Retrieved Jan 31, 2022, from https://journals.ametsoc.org/view/journals/clim/22/8/2008jcli2665.1.xml

The most recent user guide for the CMC snow depth dataset, which uses this model for the first guess field:
* https://nsidc.org/sites/nsidc.org/files/NSIDC-0447-V001-UserGuide_0.pdf
