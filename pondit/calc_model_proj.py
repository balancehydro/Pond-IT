def calc_pondit_proj(bc_calc, scalars, site, stage_storage, soils, repo_folder):
    ### calculate model using least squares best fit to calibration data
    import numpy as np
    import pandas
    

    
    params_to_fit = ['gw_out_et_percent', 'gw_out_bottom', 'rainfall_area_percent', 'soil_depth_percent', 'gw_in_percent', 'deep_fault_flow_wshed', 'deep_fault_flow_thresh', 'deep_fault_flow_lag',
                     'gw_seep_wshed', 'gw_seep_thresh', 'gw_seep_lag']
    
    
    params_model = scalars.loc[site, params_to_fit]
        

    
    t = bc_calc['mean_temp_c'].astype(np.float64)
    p = bc_calc['precip_in'].astype(np.float64)
        

    sws_calc = pondit_output(params_model, params_to_fit, scalars, p, t, bc_calc, soils, site, stage_storage)
    # sws_calc = pandas.merge(sws_calc, calib_data[['calib_wse_ft']], right_index=True, left_on='date', how='outer')

    out = params_model
    
    
    
    return sws_calc, out

def pondit_output(params_model, params_to_fit, scalars, p, t, bc_calc, soils, site, stage_storage):
    ### set up model parameters and parameters that remain constant
#     gw_out_et_percent, gw_out_bottom, rainfall_area_percent, soil_depth_percent, gw_in_percent, deep_fault_flow_wshed, deep_fault_flow_thresh, gw_seep_wshed, gw_seep_thresh, deep_fault_flow_lag, gw_seep_lag  = params_model
#     total_water_cap, wshed_area_sqft, spill_elev, spill_vol = params_cons
    import numpy as np

    scalars.loc[site, params_to_fit] = params_model

    wshed_area_sqft = scalars.loc[site, 'wshed_area_sqft']
    spill_elev = scalars.loc[site, 'spillway_elev']
    gw_out_et_percent = scalars.loc[site, 'gw_out_et_percent']
    gw_out_bottom = scalars.loc[site, 'gw_out_bottom']
    gw_in_percent = scalars.loc[site, 'gw_in_percent']
    rainfall_area_percent = scalars.loc[site, 'rainfall_area_percent']
    spill_vol = np.interp(spill_elev, stage_storage['elev'], stage_storage['storage_cuft'])
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
    sws_calc.loc[0, 'pond_elev'] = stage_storage['elev'].min()
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

            ## for months where that years total precip is larger than input threshold
            if (sws_calc.loc[n, 'precip_percent']  >= gw_seep_thresh) | (sws_calc.loc[n-1, 'gw_storage'] > 0):

                if n > np.int(gw_seep_lag * 10.0)-1:## only calculate if there are enough previous months to calculate the lag
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
                                                 stage_storage['elev'])
#         print(sws_calc)
    return sws_calc