**Contact: Aleksandra Elias Chereque // aleksandra.eliaschereque@mail.utoronto.ca**

[![DOI](https://zenodo.org/badge/454610465.svg)](https://zenodo.org/doi/10.5281/zenodo.10044950)

## Brown Temperature Index Model (B-TIM)

This simple model uses (sub-daily or daily) temperature and precipitation inputs and generates a daily record of snow depth and snow density on the same grid as the input data. Each hydrological year initializes from 0mm SWE everywhere, starting August 1 and ending July 31.

The model captures a single snow layer which evolves under the following processes:
* New snowfall -> Hedstrom/Pomeroy (1998) relationship for new snow density; snow layer density updated as weighted average
* Rain-induced melt 
* Warm-wet snow settling/densification
* Cold snow settling/densification -> follows Anderson (1978)
* Melt due to air-temperatures above -1C -> melt rate parameterization follows Kuusisto (1980)

Other parameter/constant values can be found in time_step.py.

## Instructions

* If necessary, install anaconda or miniconda (e.g., [Installation](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html))
* Create a python environment using the environment.yml file with the conda command:
```
conda env create -f environment.yml
```
* Run working example -- the last two lines printed to the console should read "TEST PASSED!"
```
conda activate env_name
bash installation_test.sh
```
* Update parameters (data_loc, latminmax, output_loc, Unique_ID) in CONFIG.py for your application
* If necessary, add details about your forcing data to utils.py (t2m_freq, tp_freq, latname, lonname, precipname, tempname)
* Adapt functions in utils.py to your forcing data and workflow.
* Verify setup by running "tests/read_files_test.ipynb". If no errors are raised, everything is ready.
* Run BTIM.py for year beginning August, YYYY and for forcing "X"
```
python BTIM.py YYYY X
```

## References \[updated Oct 2023\]
Recent publication:
* Elias Chereque et al. (2023)

Original mention of the algorithm:
* Brasnett, B. (1999). A Global Analysis of Snow Depth for Numerical Weather Prediction, Journal of Applied Meteorology, 38(6), 726-740. Retrieved Jan 31, 2022, from https://journals.ametsoc.org/view/journals/apme/38/6/1520-0450_1999_038_0726_agaosd_2.0.co_2.xml

Update:
* Ross D. Brown, Bruce Brasnett & David Robinson (2003) Gridded North American monthly snow depth and snow water equivalent for GCM evaluation, Atmosphere-Ocean, 41:1, 1-14, DOI: 10.3137/ao.410101

BTIM + snow depth observation OI used:
* Brown, R. D., & Mote, P. W. (2009). The Response of Northern Hemisphere Snow Cover to a Changing Climate, Journal of Climate, 22(8), 2124-2145. Retrieved Jan 31, 2022, from https://journals.ametsoc.org/view/journals/clim/22/8/2008jcli2665.1.xml

The most recent user guide for the CMC snow depth dataset, which uses this model for the first guess field:
* https://nsidc.org/sites/nsidc.org/files/NSIDC-0447-V001-UserGuide_0.pdf
