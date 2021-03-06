def load_data(scalars, site, folder_in):
    
    import pandas
    import datetime
    import numpy as np

   
    ## import filename referece from master input sheet
    hist_data_filename = scalars.loc[site, 'hist_data_filename'] 
    calib_area_filename = scalars.loc[site, 'calib_area_filename']
    calib_elev_filename = scalars.loc[site, 'calib_elev_filename']
    vol_area_elev_filename = scalars.loc[site, 'stage_storage_filename']
    soils_filename = scalars.loc[site, 'soils_filename']
    stage_storage_sheet_name = scalars.loc[site, 'stage_storage_sheetname']
    soils_sheet_name = scalars.loc[site, 'soils_sheetname']
    calib_area_sheet_name = scalars.loc[site, 'calib_area_sheetname']
    calib_elev_sheet_name = scalars.loc[site, 'calib_elev_sheetname']


    ### load historical climate data
    # parser = lambda dates: [datetime.datetime.strptime(x, '%Y-%m-%d') for x in dates] 
    hist_climate = pandas.read_csv(folder_in + 'hist_climate/' + hist_data_filename + '.csv', parse_dates=[0])#, date_parser=parser)
    hist_climate.columns = ['date', 'precip_in', 'mean_temp_c'] ## rename columns


    data = hist_climate.copy()
    last_hist_date = hist_climate.loc[len(hist_climate)-1, 'date']

    ## load site-specific stage-storage
    stage_storage = pandas.read_excel(folder_in + vol_area_elev_filename + '.xlsx', sheet_name=stage_storage_sheet_name, engine='openpyxl')

    
    ##load site specific soils data
    soils = pandas.read_excel(folder_in + soils_filename + '.xlsx', sheet_name=soils_sheet_name, engine='openpyxl')

    ## load calibration specified as pond area
    area_calib = pandas.read_excel(folder_in + calib_area_filename + '.xlsx', index_col='date', sheet_name=calib_area_sheet_name, engine='openpyxl')
    area_calib.index = pandas.to_datetime(area_calib.index)
    ### convert pond are to pond elevation
    area_calib['calib_wse_ft'] = np.interp(area_calib['area_sqft'], stage_storage['area_sqft'], stage_storage['elev_ft'])

    ## load elevation calibration data, recorded in feet
    elev_calib = pandas.read_excel(folder_in + calib_elev_filename + '.xlsx', index_col='date', sheet_name=calib_elev_sheet_name, engine='openpyxl')
    elev_calib.index = pandas.to_datetime(elev_calib.index)

    ## create merged calibration data 
    calib_data = pandas.merge(elev_calib, area_calib[['calib_wse_ft']], right_index=True, left_index=True, how='outer', suffixes=['_elev', '_area'])
    calib_data = calib_data.resample('MS').mean().dropna(how='all')
    calib_data['calib_wse_ft'] = calib_data.mean(skipna=True, axis=1)
    # calib_data['calib_wse_ft'] = calib_data[['calib_wse_ft_elev', 'calib_wse_ft_area']].mean(skipna=True, axis=1)

    return stage_storage, soils, data, last_hist_date, calib_data



