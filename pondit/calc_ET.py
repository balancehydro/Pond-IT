def calc_ET(data, scalars, site, folder_in):

    import scipy
    import scipy.constants
    import scipy.optimize
    import numpy as np
    import pandas
    ### calculates ET data for historical and projected time period, using Blaney Criddle calculation

    ### load blaney criddle and CIMIS data
    blaney_criddle = pandas.read_csv(folder_in + 'pondit/BlaneyCriddle_p.csv', index_col=0)
    blaney_criddle_sorted = blaney_criddle.sort_index() # sort for np.interp below
    cimis_eto_zones = pandas.read_csv(folder_in + 'pondit/CIMIS_ETO_Zones.csv', index_col=0).T
    cimis_eto_zones.index = np.int32(cimis_eto_zones.index)
    
    
    ### load lattitude and ETo zone from scalars sheet for this site
    lat = scalars.loc[site, 'lattitude']
    zone = int(scalars.loc[site, 'eto_zone'])
    
    bc_calc = data.copy()
    
    ##populate with required value types
    # bc_calc['mean_temp_c'] = scipy.constants.convert_temperature(bc_calc['mean_temp_f'], 'f', 'c')
    bc_calc['month'] = list(map(lambda x: x.month, bc_calc['date']))
    bc_calc['year'] = list(map(lambda x: x.year, bc_calc['date']))
    bc_calc['p'] = list(map(lambda x: np.interp(np.float(lat), blaney_criddle_sorted.index, 
                                                blaney_criddle_sorted.loc[:, str(bc_calc.loc[x, 'month'])]), bc_calc.index))
    bc_calc['eto_zone_' + str(zone)] = list(map(lambda x: cimis_eto_zones.loc[x, zone], bc_calc['month']))
    
    #find best values of a and b using least squares for blaney criddle calc
    def fitfunc(p, t, param): # model $f(x) = p(a*t + b)
        a,b=param
        return p *(np.multiply(a, t) + b)
    def residuals(param): # array of residuals
        return y-fitfunc(p, t, param)

    p = bc_calc['p']
    t = bc_calc['mean_temp_c'] ## use temp in C
    y = bc_calc['eto_zone_' + str(zone)] ## target ET value

    p0=[0.5, 8]#inital parameters guess
    ## fit blaney criddle parameters, using least squares
    param_best,cov,infodict,mesg,ier = scipy.optimize.leastsq(residuals, p0,full_output=True)
    print("Blaney Criddle fit parameters: a = ", param_best[0], ' b = ', param_best[1])
    

    ## calc ETo timeseries - assume ETo is good approximation of ET From pond surface using best fit values for a and b
    bc_calc['ETo'] = bc_calc['p'] * (param_best[0] * bc_calc['mean_temp_c'] + param_best[1]) ## inches/month
    
    return bc_calc