def plot_model(sws_calc, scalars, site, folder_out):

    import matplotlib.pylab as plt
    import datetime
    import numpy as np

    ## make model plot timeseries
    f, ax = plt.subplots(1, 1, figsize=(14, 6))

    sws_plot = sws_calc.loc['2000-10-01'::, :] ## assuming calibration data is not reliable prior to wy 2000 (pond sedimentation) or available
    ## plot modeled WSE
    ax.plot(np.array(sws_plot['date']), np.array(sws_plot['pond_elev']), label='Model (Historical)', zorder=2, color='C0')
 	## plot calibration data (bug in matplotlib which prevents the use of scatter with a datetime columns)
    if 'calib_wse_ft' in sws_calc.columns: #only plot calibration data for historical model which uses calibration data to predict, not for historical design runs
        ax.plot(np.array(sws_plot['date']), np.array(sws_plot['calib_wse_ft']), zorder=2, color='C1', label='Calibration Data', linestyle='', marker='o')
    

    #plot pond bottom and spillway elev
    pond_min = sws_calc.loc[:, 'pond_elev'].min().min() ## modeli initialized with pond min elev so min is always pond bottom
    pond_max = scalars.loc[site, 'spillway_elev'] 
    ax.axhline(pond_min, linestyle='--', color='k', label='Pond Bottom')
    ax.axhline(pond_max, linestyle=':', color='k', label='Pond Spillway')
    
    ##set axis limits, labels, title
    ax.set_ylim(pond_min-1, pond_max + 0.5) ## set axis limits to just include top and bottom of pond
    ax.set_xlabel('Year', fontsize=14)
    ax.set_ylabel('Pond Depth (feet)', fontsize=14)
    f.autofmt_xdate() ## render date labels diagonally
    
    ## create legend using labels
    ax.legend(loc='upper right', ncol=1, bbox_to_anchor=(1.25, 1), fontsize=14)
    ax.set_title(site, fontsize=16)
    ax.set_xlim(datetime.datetime.strptime('2000-10-01','%Y-%m-%d'), sws_calc.loc[:, 'date'].max())

    ## save figure to output folder, high quality
    f.savefig(folder_out + '/timeseries/' + site + '_timeseries_hist.png', bbox_inches='tight', dpi=330)

    ### make inflow vs outflow plot

    f1, ax1 = plt.subplots(1, 1, figsize=(14, 6))

    ## pull out flux columns
    sws_bar = sws_calc[['deep_fault_flow','gw_in','gw_in_seep','gw_out_bottom','gw_out_et', 'runoff', 'direct_rainfall', 'pond_et']].copy()
    sws_bar['gw_flux'] = sws_bar['deep_fault_flow'] + sws_bar['gw_in'] + sws_bar['gw_in_seep'] - sws_bar['gw_out_bottom'] - sws_bar['gw_out_et']
    sws_bar.index = sws_calc['date']
    sws_bar = sws_bar.loc['2005-10-01'::, :] ## assume only plot since 2005 for plotting clarity

    ## separate gw_flux positive and negative values for stacking
    sws_bar.loc[:, 'gw_flux_pos'] = sws_bar.loc[:, 'gw_flux']
    sws_bar.loc[sws_bar.loc[:, 'gw_flux']<0, 'gw_flux_pos'] = 0
    sws_bar.loc[:, 'gw_flux_neg'] = sws_bar.loc[:, 'gw_flux']
    sws_bar.loc[sws_bar.loc[:, 'gw_flux']>0, 'gw_flux_neg'] = 0

    ## plot bar chart
    w = 31 # bar width
    ### add data for stacking purposes, plot in very specific plot order
    ax1.bar(sws_bar.index, sws_bar['gw_flux_neg'] - sws_bar['pond_et'], width=w, align='center', color='xkcd:yellow orange', label='ET')
    ax1.bar(sws_bar.index, sws_bar['gw_flux_pos'] + sws_bar['direct_rainfall'] + sws_bar['runoff'], width=w, align='center', color='0.5', label='Runoff')
    ax1.bar(sws_bar.index, sws_bar['gw_flux_pos'] + sws_bar['direct_rainfall'], width=w, align='center', color='xkcd:orange', label='Direct Rainfall')
    ax1.bar(sws_bar.index, sws_bar['gw_flux_neg'], width=w, align='center', color='C0', label='Groundwater Flux')
    ax1.bar(sws_bar.index, sws_bar['gw_flux_pos'], width=w, align='center', color='C0')
    ax1.legend()

    ##set axis limits, labels, title
    ax1.set_xlabel('Year', fontsize=14)
    ax1.set_ylabel('Pond Volume (cubic feet)', fontsize=14)
    f1.autofmt_xdate() ## render date labels diagonally

    ## create legend using labels
    ax1.legend(loc='upper right', ncol=1, bbox_to_anchor=(1.25, 1), fontsize=14)
    ax1.set_title(site, fontsize=16)
    ax1.set_xlim(sws_bar.index.min(), sws_bar.index.max())

    ## save figure to output folder, high quality
    f1.savefig(folder_out + '/in_out/' + site + '_in_v_out_hist.png', bbox_inches='tight', dpi=330)

    return f, f1