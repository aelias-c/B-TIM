from os import makedirs
from os.path import exists
import sys

year = int(sys.argv[1])

forcing = str(sys.argv[2])
data_loc = 'forcing/'

mixed_pr = [0,0] #if equal, no mixed precipitation
latminmax = [40,90] 
leapdays = True

### Unique_ID will be used to name output files
Unique_ID = forcing

output_loc = 'output/'
if not exists(output_loc):
    makedirs(output_loc)
    
    
