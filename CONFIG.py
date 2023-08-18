from os import makedirs, chmod
from os.path import exists
import sys

year = int(sys.argv[1])

forcing = str(sys.argv[2])
data_loc = './input/'

mixed_pr = False
latminmax = [10,90]
lonminmax = [0,360]
leapdays = False

Unique_ID = forcing

# - For rescaling experiments - #

clim_loc = './clim/'

if len(sys.argv) > 3:
    adjust = str(sys.argv[3]) #tp, t2m, neither, both

    if adjust != 'neither':
        target_name = str(sys.argv[4])
        Unique_ID += 'r_{}_target_{}'.format(adjust, target_name)
    else:
        target_name = 'None'
else:
    target_name = forcing
    adjust = 'neither'
    
# ----------------------------- #

output_loc = './output/'+Unique_ID+'/'
if not exists(output_loc):
    makedirs(output_loc)
    
    
