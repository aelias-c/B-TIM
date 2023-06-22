from os import makedirs
from os.path import exists
import sys

year = int(sys.argv[1])

forcing = str(sys.argv[2])
data_loc = '/users/jk/20/achereque/SnowProjects2/data/01_forcing/'

mixed_pr = False
latminmax = [10,90]
lonminmax = [0,360]
leapdays = True

Unique_ID = forcing

# - For rescaling experiments - #

clim_loc = '/users/jk/20/achereque/SnowProjects2/data/01_forcing/clim/'

if len(sys.argv) > 3:
    adjust = str(sys.argv[3]) #tp, t2m, neither, both
    if adjust != 'neither':
        target_name = str(sys.argv[4])
        Unique_ID += 'r_{}_target_{}'.format(adjust, target_name)
else:
    target_name = 'None'
    adjust = 'neither'
    
# ----------------------------- #

output_loc = '/users/jk/20/achereque/SnowProjects2/data/03_reconstruction_output/'+Unique_ID+'/'
if not exists(output_loc):
    makedirs(output_loc)
    
    
