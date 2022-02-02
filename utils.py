from numpy import isin

month_names_aug = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'March', 'April', 'May', 'June', 'July']

# function calculating year-length
def len_month(m, year, leapday=True):
    '''Takes in month (Jan = 1) and year, returns the number of days in the month.'''
    if isin(m, [4, 6, 9, 11]):
        return 30
    elif m == 2:
        if (leapday) & (year % 100 == 0) & (year % 400 != 0):
            return 28
        elif (leapday) & (year % 4 == 0):
            return 29
        else:
            return 28
    else:
        return 31
    
def monthly_out_name(experiment_name, month, mixed_pr, year_tag):
        #set up save name according to settings
    savename = experiment_name + '_forced_swe_' + month + '_'

    if mixed_pr:
        savename += 'mixedpr_'
        
    savename = savename + year_tag + '.nc'
    
    return savename
