#!/usr/bin/env python3
# %%
from pathlib import Path
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec ## for figure layour

from matplotlib import rc_file

# plt.rcdefaults()
rc_file(Path(Path.home(), 'Code/covid19-vaccine/helper/matplotlibrc-custom'))

## FUNCTIONS
# %%
#### set xticks and xticklabel
def Set_Axes_Xticks(axes, start_daynum, end_daynum, data_df=None, tick_interval=30, label_fontsize=None):
    # min_daynum = start_daynum
    # max_daynum = end_daynum
    # if tick_interval < 1: tick_interval = 1
    if data_df is not None: min_daynum = min(start_daynum, min(data_df['daynum']) )
    else: min_daynum = start_daynum
    if data_df is not None: max_daynum = max(start_daynum, max(data_df['daynum']) )
    else: max_daynum = end_daynum
    if tick_interval < 1: tick_interval = 1
    
    if (max_daynum - min_daynum +1)%tick_interval == 0: xticks = np.arange(min_daynum, max_daynum+1, tick_interval) 
    else: xticks = np.arange(min_daynum, max_daynum+tick_interval, tick_interval) 
    xlabels = pd.to_datetime(xticks, unit='D', origin=pd.Timestamp('2019-12-31')).strftime('%Y-%m-%d')
    axes.set_xticks(xticks)
    if label_fontsize is not None: axes.set_xticklabels(xlabels, fontsize=label_fontsize)
    else: axes.set_xticklabels(xlabels)

    for l in axes.get_xticklabels():
        # l.set_rotation(90)
        l.set_rotation_mode('anchor')
        l.set_horizontalalignment('right')
        l.set_rotation(30)
    return axes

####
#### plot simulation trajectories and data onto a predefined axes
def Plot_On_Axes(axes, traj_arr, range_days, traj_alpha=0.03,
                tick_interval=10, iqr=True, range_95=True,
                data_df=None, data_col=None,
                title='', ylabel='' ):
    '''
    axes: matplotlib.axes
    traj_arr: numpy.ndarray with shape (nsim, nday)
    '''
    
    dat_qqqqq = np.percentile(traj_arr , q=[50, 25, 75, 2.5, 97.5], axis=0)
    # range_days = np.arange(start_daynum, end_daynum+1)


    axes.plot(range_days, (traj_arr ).T, color='grey', alpha=traj_alpha)
    axes.plot(range_days, dat_qqqqq[0,:], color='blue', alpha=1, lw=1.2)
    if iqr:
        axes.plot(range_days, dat_qqqqq[1,:]  , color='blue', alpha=1, lw=1.2)
        axes.plot(range_days, dat_qqqqq[2,:]  , color='blue', alpha=1, lw=1.2)
    if iqr and range_95:
        axes.plot(range_days, dat_qqqqq[3,:]  , color='blue', alpha=1, lw=1.2, ls='--')
        axes.plot(range_days, dat_qqqqq[4,:]  , color='blue', alpha=1, lw=1.2, ls='--')
    elif range_95 and not iqr:
        axes.plot(range_days, dat_qqqqq[3,:]  , color='blue', alpha=1, lw=1.2)
        axes.plot(range_days, dat_qqqqq[4,:]  , color='blue', alpha=1, lw=1.2)
    if data_df is not None and data_col is not None:
        axes.scatter(data_df['daynum'], data_df[ data_col ], color='black', zorder=100)
    axes = Set_Axes_Xticks(axes, range_days[0], range_days[-1], tick_interval=tick_interval, data_df=data_df)
    
    axes.set_title(title)
    axes.set_ylabel(ylabel)
    
    return axes

####
#### plot params posteriors
def Plot_All_Params_Histogram(params, id_column=True):
    offset = 1 if id_column else 0
    nparams = params.shape[1] - offset ## first column is id
    
    fig, axs = plt.subplots(ncols=4, nrows=int(nparams/4) + (nparams%4 > 0),
                            gridspec_kw={'hspace': 0.72})

    for i, ax in zip(range(nparams), axs.reshape(-1)):
        ax.hist(params.iloc[:,i+offset], bins=40)
        ax.set_title(params.columns[i+offset], fontsize=15, pad=6)
    
    return fig, axs

#### plot params posteriors
def Plot_All_Params_Samples(params, id_column=True, num_samples=-1):
    offset = 1 if id_column else 0
    nparams = params.shape[1] - offset ## first column is id
    if num_samples < 1: num_samples=params.shape[0]
    ncols =2
    nrows = int(nparams/ncols) + (nparams%ncols > 0)
    
    # fig, axs = plt.subplots(ncols=1, nrows=nparams, figsize=(15, 4*nparams))
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows,
                            figsize=(15*ncols, 4*nrows))

    for i, ax in zip(range(nparams), axs.reshape(-1)):
        ax.plot((params.iloc[-num_samples:,i+offset]).T, color='grey', marker='x', ls='--', alpha=min(1.0, 100.0/num_samples))
        ax.set_title(params.columns[i+offset], fontsize=15, pad=6)
    
    return fig, axs

####
#### plot betas posterior, every 7 days
#### these are equivalent beta spline coefficients
def Plot_Weekly_Betas_Histogram(betas, id_column=True):
    offset = 1 if id_column else 0
    nbetas = betas.shape[1] - offset ## first column is id
    
    nweeks = int(nbetas/7) + 1 ## +1 for the first beta

    fig, axs = plt.subplots(ncols=4, nrows=int(nweeks/4) + (nweeks%4 > 0), sharex='row',
                            gridspec_kw={'hspace': 0.72})

    for i, ax in zip(range(0, nbetas+1, 7), axs.reshape(-1)):
        ax.hist(betas.iloc[:,i+offset], bins=40)
        ax.set_title('beta {}'.format(betas.columns[i+offset]), fontsize=15, pad=6)
    
    return fig, axs

####
#### plot all betas, with median line, IQR and 95% credible range
def Plot_All_Betas(betas_posterior, start_daynum, end_daynum, tick_interval=10, iqr=True, range_95=True, days_preepidemics=10,
                   title='', ylabel='' ):
    if tick_interval < 1: tick_interval = 1
    # if days_preepidemics < 1: days_preepidemics = 1
    if days_preepidemics > betas_posterior.shape[1]: days_preepidemics = betas_posterior.shape[1]

    range_days = range(start_daynum, end_daynum+1)
    if days_preepidemics < 1: preepi_mean = 1.0
    # else: preepi_mean = np.mean( np.mean(betas_posterior.iloc[:, :days_preepidemics], axis=1) )
    else: preepi_mean = betas_posterior.iloc[:, :days_preepidemics].mean().mean()
    betas_qqqqq = betas_posterior.quantile(q=[0.5, 0.25, 0.75, 0.025, 0.975])

    fig, ax = plt.subplots()
    ax = Set_Axes_Xticks(ax, start_daynum, end_daynum, tick_interval=tick_interval)

    ax.plot(range_days, (betas_posterior / preepi_mean).T, color='grey', alpha=0.04)
    ax.plot(range_days, betas_qqqqq.iloc[0,:] / preepi_mean , color='blue', alpha=1, lw=1.2)
    if iqr:
        ax.plot(range_days, betas_qqqqq.iloc[1,:] / preepi_mean , color='blue', alpha=1, lw=1.2)
        ax.plot(range_days, betas_qqqqq.iloc[2,:] / preepi_mean , color='blue', alpha=1, lw=1.2)
    if iqr and range_95:
        ax.plot(range_days, betas_qqqqq.iloc[3,:] / preepi_mean , color='blue', alpha=1, lw=1.2, ls='--')
        ax.plot(range_days, betas_qqqqq.iloc[4,:] / preepi_mean , color='blue', alpha=1, lw=1.2, ls='--')
    elif range_95 and not iqr:
        ax.plot(range_days, betas_qqqqq.iloc[3,:] / preepi_mean , color='blue', alpha=1, lw=1.2)
        ax.plot(range_days, betas_qqqqq.iloc[4,:] / preepi_mean , color='blue', alpha=1, lw=1.2)

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    
    return fig, ax

# %%
if __name__ == "__main__":
    pass

