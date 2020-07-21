def load_climate_data(data_in, model, scalars, site, folder_in, proj_pre_raw, proj_temp_raw, last_hist_date):
    
    import pandas
    import numpy as np
    import calendar
    climate_data = pandas.DataFrame(proj_pre_raw['date'])
    climate_data['num_days_month'] = list(map(lambda x: np.float(calendar.monthrange(x.year, x.month)[1]), climate_data['date']))
    climate_data['precip_in_proj'] = proj_pre_raw.loc[:, model] * climate_data['num_days_month']  / 25.4 ## convert from mm/day to inches/month
    climate_data['mean_temp_f_proj'] = (proj_temp_raw.loc[:, model] * (9.0/5.0)) + 32.0 ### convert from C to F

    
    data_merged = data_in.merge(climate_data[['date', 'precip_in_proj', 'mean_temp_f_proj']], on='date', how='outer') ## join climate data and historical data 
    data_merged['month'] = list(map(lambda x: x.month, data_merged['date']))
    data_merged['year'] = list(map(lambda x: x.year, data_merged['date']))
    data_merged['wy'] = data_merged['year']
    data_merged.loc[data_merged['month'].isin([10, 11, 12]), 'wy'] = data_merged.loc[data_merged['month'].isin([10, 11, 12]), 'year'] + 1
    data_merged['count_col'] = data_merged['mean_temp_f']

    ### limit to complete water years
    annual_climate_hist = data_merged.groupby('wy').agg({'precip_in': 'sum', 'precip_in_proj':'sum', 'mean_temp_f':'mean', 'mean_temp_f_proj':'mean', 'count_col':'count'}).loc[1976:2017]
    annual_climate_hist = annual_climate_hist.loc[annual_climate_hist['count_col'] ==12, :]

    precip_scale = annual_climate_hist['precip_in'].median() / annual_climate_hist['precip_in_proj'].median()
    temp_scale = annual_climate_hist['mean_temp_f'].median() / annual_climate_hist['mean_temp_f_proj'].median()
    
    data_merged['precip_in_proj'] = data_merged['precip_in_proj'] * precip_scale
    data_merged['mean_temp_f_proj'] = data_merged['mean_temp_f_proj'] * temp_scale

    data_out = data_merged.copy()
    data_out['precip_in'].fillna(data_out['precip_in_proj'], inplace=True)
    data_out['mean_temp_f'].fillna(data_out['mean_temp_f_proj'], inplace=True)
    data_out.drop(['mean_temp_f_proj', 'precip_in_proj'], axis=1, inplace=True)
    
    return data_out
    
def load_all_climate_data(scalars, site, folder_in):
    
    import pandas
    pre_proj_filename = scalars.loc[site, 'proj_precip_filename']
    temp_proj_filename = scalars.loc[site, 'proj_temp_filename']
    
    
    proj_pre_raw = pandas.read_csv(folder_in + 'proj_climate/' + pre_proj_filename + '.csv', parse_dates=['date'])
    proj_temp_raw = pandas.read_csv(folder_in + 'proj_climate/' + temp_proj_filename + '.csv', parse_dates=['date'])
    return proj_pre_raw, proj_temp_raw