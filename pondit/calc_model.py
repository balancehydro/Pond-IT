def calc_pondit(bc_calc, scalars, site, stage_storage, soils, repo_folder, calib_data):
    import pandas
    import numpy as np
    import scipy
    import scipy.optimize

    ### load model parameters from scalars sheet for model calculations
    default_params = pandas.read_csv(repo_folder + 'pondit/' + 'Default_Parameters.csv', index_col=0).T
    default_params.index = [x if x!='value' else site for x in default_params.index]

    ## divide lag times by 10 so they are in better numerical scale for model prediction
    default_params.loc[:, ['deep_fault_flow_lag', 'gw_seep_lag']]  = default_params.loc[:, ['deep_fault_flow_lag', 'gw_seep_lag']] / 10.0 

    scalars = pandas.merge(scalars, default_params, right_index=True, left_index=True) ## master inputs with default params for predicting model
    
    ## load indicators for running either gw modules
    gw_seep = scalars.loc[site, 'gw_seep']
    gw_fault = scalars.loc[site, 'gw_fault']
    
    ## complete list of potential model fit parameters
    params_to_fit = list(default_params.columns)

    ## remove gw_seep parameters if they shouldn't be included in predicting the model
    if gw_seep == 'n':
        gw_seep_params = ['gw_seep_wshed', 'gw_seep_thresh', 'gw_seep_lag']
        params_to_fit = [ x for x in params_to_fit if x not in (gw_seep_params)]
        scalars.loc[site, gw_seep_params] = 0 

    ## remove gw_fault parameters if they shouldn't be included in predicting the model
    if gw_fault == 'n':
        gw_fault_params = ['deep_fault_flow_wshed', 'deep_fault_flow_thresh', 'deep_fault_flow_lag']
        params_to_fit = [ x for x in params_to_fit if x not in (gw_fault_params)]
        scalars.loc[site, gw_fault_params] = 0 
    
    ## parameters to not fit (ie constants) are all parameters not including params to fit
    params_not_fit = [x for x in list(default_params.columns) if x not in params_to_fit]
    

    bc_calc_all = bc_calc.copy() ## separate all data from model fit period
    bc_calc_model = bc_calc.loc[bc_calc['date'] >= '2000-10-01', :] ## assumed unreliable or unavailable calibration prior to 2000
    bc_calc_model.index = np.arange(0, len(bc_calc_model), 1)

    params_model = scalars.loc[site, params_to_fit]
    params_cons = scalars.loc[site, params_not_fit]

    ## define function to fit model
    def fitfunc(params_model, params_to_fit, scalars, bc_calc_model, soils): # model 
        y = pondit_calib(params_model, params_to_fit, scalars, bc_calc_model, soils, calib_data, site, stage_storage)
        return y
    ## define function to calculate residuals to minimize
    def residuals(params_model): 
        import numpy as np
        return np.power(y-fitfunc(params_model, params_to_fit, scalars, bc_calc_model, soils), 2).sum()
        
    ## calibration data
    y = np.array(calib_data['calib_wse_ft'])

    p0=params_model #inital parameters guess

    ## define bounds using default parameters sheet
    bounds = ()
    for param in params_to_fit:
        bounds = bounds + ((default_params.loc['min', param], default_params.loc['max', param]),)
    
    ## minimize residuals
    out = scipy.optimize.minimize(residuals, p0, method='SLSQP', bounds=bounds, options={'disp':True, 'eps':0.1, 'ftol':0.005})
    print(out) ## print final result

    ## calculate final results with optimized model parameters for entier historical period for output
    sws_calc = pondit_output(out.x, params_to_fit, scalars, bc_calc_all, soils, site, stage_storage)
    sws_calc = pandas.merge(sws_calc, calib_data[['calib_wse_ft']], right_index=True, left_on='date', how='outer') ## add calibration data to output
    
    ## save optimal model parameters
    scalars_out = scalars.copy()
    scalars_out.loc[site, params_to_fit] = out.x
    scalars_out.loc[:, ['deep_fault_flow_lag', 'gw_seep_lag']]  = np.int64(scalars_out.loc[:, ['deep_fault_flow_lag', 'gw_seep_lag']] * 10.0) # multiple lag by 10 to put into month units
    
    
    return sws_calc, out, scalars_out


def pondit_calib(params_model, params_to_fit, scalars, bc_calc_model, soils, calib_data, site, stage_storage):
    
    import numpy as np

    sws_calc = pondit_output(params_model, params_to_fit, scalars, bc_calc_model, soils, site, stage_storage)
    
    
    y = calib_data['calib_wse_ft']

    out = np.array(sws_calc.loc[sws_calc['date'].isin(y.index), 'pond_elev'])
    print(params_model, 'sum squared error: ', np.power(np.array(y)-out, 2).sum())
    return out

def pondit_output(params_model, params_to_fit, scalars, bc_calc, soils, site, stage_storage):
    ### set up model parameters and parameters that remain constant
    import numpy as np

    scalars.loc[site, params_to_fit] = params_model

    wshed_area_sqft = scalars.loc[site, 'wshed_area_sqft']
    spill_elev = scalars.loc[site, 'spillway_elev']
    gw_out_et_percent = scalars.loc[site, 'gw_out_et_percent']
    gw_out_bottom = scalars.loc[site, 'gw_out_bottom']
    gw_in_percent = scalars.loc[site, 'gw_in_percent']
    rainfall_area_percent = scalars.loc[site, 'rainfall_area_percent']
    spill_vol = np.interp(spill_elev, stage_storage['elev_ft'], stage_storage['storage_cuft'])
    zone = int(scalars.loc[site, 'eto_zone'])
    gw_seep_lag = scalars.loc[site, 'gw_seep_lag']
    gw_seep_thresh = scalars.loc[site, 'gw_seep_thresh']
    gw_seep_wshed = scalars.loc[site, 'gw_seep_wshed']
    deep_fault_flow_wshed = scalars.loc[site, 'deep_fault_flow_wshed']
    deep_fault_flow_lag = scalars.loc[site, 'deep_fault_flow_lag']
    deep_fault_flow_thresh = scalars.loc[site, 'deep_fault_flow_thresh']
    soil_depth_percent = scalars.loc[site, 'soil_depth_percent']
    
    
    ## drop un-used columns for simplicity
    sws_calc = bc_calc.copy()

    ## calculate various soil metrics: total water capacity, total root zone capacity
    percent_impervious_cov = (1 - soils['Percent of AOI'].sum() )
    area_weight_water_cap = ((soils['Percent of AOI'] * soils['Depth Weighted Water Capacity']).sum() /
                                soils['Percent of AOI'].sum()) #in/in
    area_weight_profile_depth = ((soils['Percent of AOI'] * soils['Profile Thickness (inches)']).sum() / 
                                     (1 - percent_impervious_cov)) #inches
    total_water_cap = area_weight_profile_depth * area_weight_water_cap  * soil_depth_percent #inches
    root_water_cap = area_weight_water_cap * 18.0 ## assume 18-inch root zone


    ## calculate water year from year, start with copying calendar year
    sws_calc['wy'] = sws_calc['year']
    ## if month is oct, nov, dec, add one to calendar year to get water year
    sws_calc.loc[sws_calc['month'].isin([10, 11, 12]), 'wy'] = np.int64(sws_calc.loc[sws_calc['month'].isin([10, 11, 12]), 'year']) +1

    ## Calculate metrics that aren't iterative
    sws_calc['annual_precip'] = list(map(lambda x: sws_calc.loc[sws_calc['wy'] == x, 'precip_in'].sum(), sws_calc['wy'])) ## total precip by wy
    sws_calc['mean_annual_precip'] = sws_calc.loc[sws_calc['month']==1, 'annual_precip'].mean() ## average historical precip, all wy (e.g. 'normal')
    sws_calc['precip_percent'] = sws_calc['annual_precip'] / sws_calc['mean_annual_precip'] ## wy precip / annual average precip
    ## Calculate last years precip amount ratio relative to 'normal'
    sws_calc['precip_percent_lag1'] = 1 ## fill in default value of 1 for first year of record, serves as memory effect, higher when last year was wet, etc
    sws_calc.loc[sws_calc['wy'] > sws_calc['wy'].min(), 'precip_percent_lag1'] = sws_calc['precip_percent'].shift(12) ## lag one year to use last years ratio in calcs

    ## calc cumulative soil water
    sws_calc.loc[0, 'soil_water'] = 0 ## initialize timeseries
    sws_calc.sort_index(inplace=True)
    
    for n in sws_calc.index[1::]:
        ## soil water = last months soil water + precip - ETo, with a max of total water capacity and min of zero
        sws_calc.loc[n, 'soil_water'] = max(min(sws_calc.loc[n-1, 'soil_water'] + sws_calc.loc[n, 'precip_in'] - ### was p[n] here
                                                sws_calc.loc[n, 'ETo'], total_water_cap), 0)

    ## calc how much was went into the soil that month
    sws_calc['water_into_soil'] = sws_calc['soil_water'].diff(1).clip(lower=0)
    sws_calc.loc[0, 'water_into_soil'] = 0 ##initialize
    ## calculate excess water for that month: precip - water that soakied into soil - ETo, min of zero
    sws_calc['excess_water'] = (sws_calc['precip_in'] - sws_calc['water_into_soil'] - sws_calc['ETo']).clip(lower=0) #inches

    ## calc runoff: excess water(in feet) times watershed area
    sws_calc['runoff'] = ((sws_calc['excess_water']) / 12.) * wshed_area_sqft #cubic feet

    ## calculate how full soil column is, percentage
    sws_calc['water_cap_percent'] = (1.0 - (sws_calc['soil_water'] / total_water_cap))



    ## calculate deep fault flow: rolling average of precip (over 6 months), shifted by specified lag, scaled by input parameter, over the whole watershed
    sws_calc['deep_fault_flow'] = 0
    if np.int(deep_fault_flow_lag * 10.0) > 0:
        sws_calc['deep_fault_flow'] = ((sws_calc['precip_in'].rolling(window=6,center=True).mean().shift(int(deep_fault_flow_lag * 10.0))) 
                                       / 12.0 * wshed_area_sqft * deep_fault_flow_wshed).fillna(0)
        
        if deep_fault_flow_thresh != 0:
            sws_calc.loc[sws_calc['precip_percent'] < deep_fault_flow_thresh, 'deep_fault_flow'] = 0
        


    ## initialize various parameters
    sws_calc.loc[0, 'pond_et'] = 0
    sws_calc.loc[0, 'pond_area'] = 0
    sws_calc.loc[0, 'pond_volume'] = 0
    sws_calc.loc[0, 'pond_elev'] = stage_storage['elev_ft'].min()
    sws_calc.loc[0, 'direct_rainfall'] = 0
#     sws_calc.loc[0, 'not_root_water'] = 0
    sws_calc['gw_storage'] = 0
    sws_calc['gw_in_seep'] = 0
    sws_calc['gw_seep_percent'] = 0
    sws_calc.loc[0, 'gw_in'] = 0
    sws_calc['gw_out_et'] = 0
    sws_calc['gw_out_bottom'] = 0

    ## calculate maximum gw storage of the watershed below the root zone
    max_gw_storage = (total_water_cap - root_water_cap) / 12.0 * (gw_seep_wshed * wshed_area_sqft) #cu ft
    ## max pond volume
    max_pond_vol = stage_storage['storage_cuft'].max()

    ## calculate how much water is in soil column, but below root zone
    sws_calc['not_root_water'] = (sws_calc['soil_water'] - sws_calc['ETo'] - root_water_cap).clip(lower=0) / 12.0 * gw_seep_wshed * wshed_area_sqft #cu ft


    ## loop through whole time series (historical and projected)
    for n in sws_calc.index[1::]:
        ## calculate direct rainfall: last months pond area plus pond fringe area times monthly precip in feet
        sws_calc.loc[n, 'direct_rainfall'] = ((sws_calc.loc[n-1, 'pond_area'] * rainfall_area_percent) * 
                                              (sws_calc.loc[n, 'precip_in'] / 12.)) ## cu ft


        ## calculate groundwater output as a function of ET - acts to adjust ET loses from aquatic veg, tree cover, soil properties, soil moisutre wicking.
        ## last months ET, over the pond fringe area, scaled by how saturated the soil column is, and the model input parameter
        sws_calc.loc[n, 'gw_out_et'] = (sws_calc.loc[n-1, 'pond_et'] * rainfall_area_percent * gw_out_et_percent * sws_calc.loc[n, 'water_cap_percent'])

        ## calculate groundwater losses because pond is 'leaky' - leaky berm, high permeability soils, rapid pathways to fractures, etc
        ## function of pond volume (pressure drain more volume from pond): last months pond volume times model input parameter
        sws_calc.loc[n, 'gw_out_bottom'] = (sws_calc.loc[n-1, 'pond_volume'] * gw_out_bottom)
        
        ### if model parameterization has groundwater input seep: the seep module represents groundwater inputs from shallow bedrock fractures, 
        ### or other medium-term lagged groundwater inputs that may contribute to seeps
        if np.int(gw_seep_lag * 10.0) > 0:

            ## for months where that year's total precip is larger than input threshold
            if (sws_calc.loc[n, 'precip_percent']  >= gw_seep_thresh) | (sws_calc.loc[n-1, 'gw_storage'] > 0):

                if n > np.int(gw_seep_lag *10.0)-1:## only calculate if there are enough previous months to calculate the lag
                    ## calculate gw storage below root zone: last months gw storage + the lagged 'not root water' (i.e. water stored below root zone)
                    sws_calc.loc[n, 'gw_storage'] = sws_calc.loc[n-1, 'gw_storage'] + sws_calc.loc[n - np.int(gw_seep_lag * 10.0), 'not_root_water']

                    ## calculate the percent by which to discharge the stored groundwater into the pond, 
                    ## calculated as percent full, rescaled to run from 25 to 92%, (so that all of it can't leave at once, and so something always discharges if available)
                    sws_calc.loc[n, 'gw_seep_percent'] = (sws_calc.loc[n, 'gw_storage'] / max_gw_storage) / 1.5 + 0.25

                    ## calculate min discharge from seep: min of either last months gw storage or last months losses (both gw out and ET)
                    ## i.e. can't discharge more water than leaves the pond to control the rate of groundwater discharge
                    gw_in_seep_min = min(sws_calc.loc[n-1, 'gw_storage'], sws_calc.loc[n-1, 'pond_et'] + sws_calc.loc[n-1, 'gw_out_et'] + sws_calc.loc[n-1, 'gw_out_bottom'])

                    ## calculate seep discharge amount: gw storage scaled by input parameter percent, min set to above parameter
                    sws_calc.loc[n, 'gw_in_seep'] = (sws_calc.loc[n, 'gw_storage'] * sws_calc.loc[n, 'gw_seep_percent']).clip(gw_in_seep_min)

                    ## update this months gw storage (reduce by the amount discharged)
                    sws_calc.loc[n, 'gw_storage'] = (sws_calc.loc[n, 'gw_storage'] - sws_calc.loc[n, 'gw_in_seep']).clip(0)


        ## calculate groundwater input via slow infiltration from the pond fringe area: last months direct rainfall scaled by model input parameter and the 
        sws_calc.loc[n, 'gw_in'] = sws_calc.loc[n-1, 'direct_rainfall'] * gw_in_percent * sws_calc.loc[n, 'precip_percent_lag1']

        ## calculate all inputs
        sws_calc.loc[n, 'ins'] = (sws_calc.loc[n, 'runoff'] + sws_calc.loc[n, 'gw_in'] + sws_calc.loc[n, 'direct_rainfall'] + 
                                  sws_calc.loc[n, 'deep_fault_flow'] + sws_calc.loc[n, 'gw_in_seep'])
        ## calculate all outputs: use last months ET to calculate pond volume
        ## this is because pond area is needed to calculate ET, which need pond volume, so either use last months pond area, or last months et; chose the latter
        sws_calc.loc[n, 'outs'] = sws_calc.loc[n, 'gw_out_bottom'] + sws_calc.loc[n-1, 'pond_et'] + sws_calc.loc[n, 'gw_out_et']

        ##update pond volume: previous months volume + inputs - outputs
        sws_calc.loc[n, 'pond_volume'] = (sws_calc.loc[n-1, 'pond_volume'] +  sws_calc.loc[n, 'ins'] - sws_calc.loc[n, 'outs']).clip(0, spill_vol)

        ## use stage-storage to infer pond area
        sws_calc.loc[n, 'pond_area'] = np.interp(sws_calc.loc[n, 'pond_volume'], stage_storage['storage_cuft'], 
                                                 stage_storage['area_sqft'])
        ## calculate this months ET with updated pond area, interpolated from stage-storage
        sws_calc.loc[n, 'pond_et'] = sws_calc.loc[n, 'pond_area'] * (sws_calc.loc[n, 'ETo'] / 12.)

        #calculate this months pond elevation, interpolated from stage-storage
        sws_calc.loc[n, 'pond_elev'] = np.interp(sws_calc.loc[n, 'pond_volume'], stage_storage['storage_cuft'], 
                                                 stage_storage['elev_ft'])
#         print(sws_calc)
    return sws_calc