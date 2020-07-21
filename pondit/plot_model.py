def plot_model(sws_calc, scalars, site, folder_out):

    import matplotlib.pylab as plt
    

    ## make model plot timeseries
    f, ax = plt.subplots(1, 1, figsize=(14, 6))

    sws_plot = sws_calc.loc['2000-10-01'::, :] ## assuming calibration data is not reliable prior to wy 2000 (pond sedimentation) or available
    ## plot modeled WSE
    ax.plot(sws_plot['date'], sws_plot['pond_elev'], label='Model (Historical)', zorder=2, color='C0')
 	## plot calibration data (bug in matplotlib which prevents the use of scatter with a datetime columns)
    ax.plot(sws_plot['date'], sws_plot['calib_wse_ft'], zorder=2, color='C1', label='Calibration Data', linestyle='', marker='o')
    

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
    ax.set_xlim('2000-10-01', sws_calc.loc[:, 'date'].max())

    ## save figure to output folder, high quality
    f.savefig(folder_out + '/timeseries/' + site + '_timeseries_hist.png', bbox_inches='tight', dpi=330)

    return f