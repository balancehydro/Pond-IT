

def plot_hydroperiod(sws_calc, scalars, site, folder_out):
    import pandas
    import numpy as np
    import matplotlib.pylab as plt
    import glob
    from matplotlib.ticker import FuncFormatter
    import matplotlib.colors
    import seaborn

    ## pond minimum pond elevation to indicate when pond is dry
    pond_min = sws_calc.loc[:, :, 'pond_elev'].min().min() ## model calculation always initialized with minimum pond elevation

    plot_inundation = scalars.loc[site, 'plot_inundation'] ## target ponding depth for results plot
    
    ### define colormap
    colors_blue = ['#e9f4fb', '#d3e9f8', '#a8d2f0', '#7cbce9', '#51a5e1', '#258fda', '#1e72ae', '#165683', '#0f3957', 'black']
    cm_blue = matplotlib.colors.LinearSegmentedColormap.from_list(
            'Blues', colors_blue, N=10)
    

    ## pull out pond elevation data
    pond_elev = sws_calc.loc[:, :, 'pond_elev'].copy()
    pond_elev.index = sws_calc.loc[sws_calc.items[0],  :, 'date'] ## specify index as date

    ## target depth is pond_min plus plot_inundation
    pond_elev[pond_elev <= pond_min + plot_inundation] = -9999 ## when depth less than target, set to nan
    pond_elev[pond_elev > pond_min + plot_inundation] = 1  ### when greater than target, set to 1 ("inundated")
    pond_elev[pond_elev == -9999] = 0 ### set nan values to zero ("not inundated") - do this so don't unintentially override actual elevation data

    ## Calculate percent inundation 
    pond_percent = pandas.DataFrame(np.round((pond_elev).mean(1)* 10) / 10, index=pond_elev.index, columns=['percent_inundation'])
    ## calculate water year
    pond_percent['month'] = list(map(lambda x: x.month, pond_percent.index)) 
    pond_percent['year'] = list(map(lambda x: x.year, pond_percent.index))
    pond_percent['wy'] = pond_percent['year']
    pond_percent.loc[pond_percent['month'].isin([10, 11, 12]), 'wy'] = pond_percent.loc[pond_percent['month'].isin([10, 11, 12]), 'year'] + 1
    ## calculate decade start
    pond_percent['decade_start'] = np.int_(np.floor(pond_percent['wy'] / 10) * 10)

    ##pivot data for to arrange in grid (month x wy)
    pivot = pandas.pivot_table(pond_percent, values='percent_inundation', index='wy', columns='month')
    pivot = pivot[[10, 11, 12, 1, 2, 3, 4, 5, 6, 7,8, 9]] ## reorder columns to be in wy sequence
    pivot.index.name = 'Water Year' ## rename column so it prints on heatmap
    
    ### create water year grid inundation plot
    f, ax = plt.subplots(1, 1, figsize=(4, 10)) 
    pivot.columns.name = 'Month'
    seaborn.heatmap(pivot, cmap=cm_blue, ax=ax, cbar=True, linewidths=1, vmin=0, vmax=1, annot_kws={"fontsize": 60})
    ax.set_title(str(plot_inundation) + '-foot Inundation Probability', fontsize=12)
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=8)
    plt.tight_layout()
    f.savefig(folder_out + '/grid plots/'+ site + '.png', bbox_inches='tight', dpi=330) # save to results folder
    
    ## pivot data again but by decade
    decade_pivot = pandas.pivot_table(pond_percent, values='percent_inundation', index='decade_start', columns='month', aggfunc='mean')
    decade_pivot = decade_pivot[[10, 11, 12, 1, 2, 3, 4, 5, 6, 7,8, 9]] ## reorder columns to be in wy sequence
    decade_pivot.index = list(map(lambda x: str(x) + ' - ' + str(x + 10), decade_pivot.index))

    ## Create decade plot
    f, ax = plt.subplots(1, 1, figsize=(6, 5)) 
    decade_pivot.columns.name = 'Month'
    fmt = lambda x,pos: '{:.0%}'.format(x) ## format annotations
    seaborn.heatmap(decade_pivot, cmap=cm_blue, ax=ax, cbar=True, linewidths=1, vmin=0, vmax=1, annot=True, fmt='.0%', cbar_kws={'format': FuncFormatter(fmt)}, annot_kws={"fontsize":6})
    ax.set_title(str(plot_inundation) + '-foot Inundation Probability', fontsize=12)
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=8)

    ## save to results folder
    f.savefig(folder_out + '/grid plots decade/'+ site + '_' + '.png', bbox_inches='tight', dpi=300)
    
    ## turn off seaborn so it doesn't affect other plots
    seaborn.reset_orig()
    plt.close('all')
    
    return f
