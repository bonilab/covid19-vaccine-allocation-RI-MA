#!/usr/bin/env python3
# %%
import numpy as np
# import pandas as pd 
from pathlib import Path
from cycler import cycler
import matplotlib.pyplot as plt 
from matplotlib import cm
# from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import rc_file

# plt.rcdefaults()
rc_file(Path(Path.home(), 'Code/covid19-vaccine/helper/matplotlibrc-custom'))

# %%
def Get_Efficacy(time, half_life, slope):
    return (half_life**slope) / ((half_life**slope) + (time**slope)) 

# %%
if __name__ == "__main__":
    # %%
    base_dir = Path(Path.home(), 'Code/covid19-vaccine/helper/')
    fz=25 ## font size

    # %%
    half_life = 180 ## in day 
    slope = 2
    ndays = 365
    foi = 0.3
    eff = [Get_Efficacy(x, half_life, slope) for x in range(ndays)]
    ckpts = [x*30 for x in (2,3,4,5,6)]
    if half_life not in ckpts and half_life < ndays: 
        ckpts = ckpts + [int(half_life)]

    
    # %%
    ####
    #### vaccine profile screening
    ####
    hl_ar = [180, 360, 540]
    len_hl_ar = len(hl_ar)
    sl_ar = [1.5, 2.0, 3.0, 4.0]
    # sl_ar = [2.0, 3.0, 4.0] # 180
    # sl_ar = [1.5, 2.0] # 360
    # sl_ar = [1.5, 2.0] # 540
    len_sl_ar = len(sl_ar)
    ndays = 1096
    eff_ar = np.array(list(map(lambda hl, sl: [Get_Efficacy(x, hl, sl) for x in range(ndays)], 
                            np.repeat(hl_ar, len_sl_ar), 
                            np.tile(sl_ar, len_hl_ar) )))

    # %%
    #### plot and save separate profiles
    ndays = 391
    for half_life in hl_ar:
        for slope in sl_ar:
            ckpts = [x*30 for x in (2,3,4,5,6,12)]
            eff = [Get_Efficacy(x, half_life, slope) for x in range(ndays)]
            # if half_life not in ckpts and half_life < ndays: 
            #     ckpts = ckpts + [int(half_life)]

            #### plot one single profile
            fig, ax = plt.subplots(figsize=(14,9))
            ax.plot(range(ndays), eff)
            # ax.axvline(x=half_life, ls='--', lw=0.8, color='grey')
            for t in ckpts:
                ax.annotate('D{}\n{:.1%}'.format(t,eff[t]), xy=(t, eff[t]), 
                            xytext=(t, eff[t] -0.05),
                            ha="center", va='top', fontsize=18, 
                            arrowprops=dict(arrowstyle='simple',
                                            facecolor='tab:red',
                                            color='tab:red'))
            # ax.text(ndays, 0.98, 'half-life: {}\nslope: {}'.format(half_life, slope), 
            #         ha='right', va='top', wrap=True, fontsize=11,
            #         bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
            ax.margins(y=0.12)
            ax.set_title('vaccine half-life: {}; vaccine slope: {}'.format(half_life, slope))
            # ax.set_xticks([x for x in range(0, ndays+1, 30)] + [half_life] )
            # ax.set_xlim(0,365)
            ax.set_xticks(range(0, ndays+1, 30))
            ax.set_xlabel('Days since vaccination')
            ax.set_yticks(np.arange(0, 1.01, 0.2))
            ax.set_yticklabels(['{:.0%}'.format(x) for x in  np.arange(0, 1.01, 0.2)])
            ax.set_ylabel('Vaccince efficacy')
            
            fig.savefig(Path(base_dir, 'efficacy_curve-halflife_{}_slope_{}.png'.format(half_life, slope)))


    # %%
    cl_tab10 = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
                'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    ls_mix = ['-', '-.', ':', '--', (0,(1,5,1,5,3,5)), (0,(1,5,1,5,3,5,3,5)), (0, (1, 10)) ]
    # style_cycler = cycler(color=cl_tab10[:len_sl_ar], linestyle=ls_mix[:len_hl_ar])
    style_cycler = cycler(color=cl_tab10[:len_sl_ar])

    plt.rc('axes', prop_cycle=style_cycler )
    # plt.rcdefaults()

    # %%
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20,10), dpi=120,
                           gridspec_kw={'width_ratios': [15, 1], 'wspace': 0.05, 'right': 0.87})
    # sl_ar = [2.0, 3.0, 4.0] # 180
    # sl_ar = [1.5, 2.0] # 360
    # sl_ar = [1.5, 2.0] # 540
    for i in range(len_hl_ar):
        # axs[0].plot(eff_ar[(i*len_sl_ar):((i+1)*len_sl_ar ),:].T, ls=ls_mix[i])
        if hl_ar[i] == 180:
            axs[0].plot(eff_ar[(i*len_sl_ar),:], ls=ls_mix[i], alpha=0.0)
            axs[0].plot(eff_ar[(i*len_sl_ar+1):((i+1)*len_sl_ar ),:].T, ls=ls_mix[i])
        elif hl_ar[i] == 360:
            axs[0].plot(eff_ar[(i*len_sl_ar):((i+1)*len_sl_ar -2),:].T, ls=ls_mix[i])
            axs[0].plot(eff_ar[((i+1)*len_sl_ar -2):((i+1)*len_sl_ar),:].T, ls=ls_mix[i], alpha=0.0)
        elif hl_ar[i] == 540:
            axs[0].plot(eff_ar[(i*len_sl_ar):((i+1)*len_sl_ar -2),:].T, ls=ls_mix[i])
            # axs[0].plot(eff_ar[((i+1)*len_sl_ar -2):((i+1)*len_sl_ar),:].T, ls=ls_mix[i], alpha=0.0)


    axs[0].axhline(y=0.90, color='grey')
    axs[0].axvline(x=60, color='grey')
    axs[0].axvline(x=90, color='grey')
    # axs[0].axhline(y=0.95, color='grey')
    # axs[0].axvline(x=40, color='grey')
    # axs[0].axvline(x=120, color='grey')
    axs[0].set_yticks(np.arange(0, 1.01, 0.1))
    axs[0].set_yticklabels(['{:.0%}'.format(x) for x in  np.arange(0, 1.01, 0.1)])
    # axs[0].set_xticks(np.arange(0, 181, 30))
    # axs[0].set_xlim(0,180)
    # axs[0].set_ylim(0.80, 1.0)
    axs[0].set_xlabel('Days since vaccination')
    axs[0].set_ylabel('Vaccince efficacy')
    axs[0].margins(x=0.03, y=0.03)

    for i in range(len_hl_ar):
        axs[1].plot([0,1], [len_sl_ar+len_hl_ar-i, len_sl_ar+len_hl_ar-i], lw=4, ls=ls_mix[i], color='black')
        axs[1].text(1.2, len_sl_ar+len_hl_ar-i, 'half-life {}'.format(hl_ar[i]), 
                    fontsize=fz, ha='left', va='center')
    for i in range(len_sl_ar):
        axs[1].plot([0,1], [len_sl_ar-i, len_sl_ar-i], lw=4, ls='-')
        axs[1].text(1.2, len_sl_ar-i, 'slope {}'.format(sl_ar[i]), 
                    fontsize=fz, ha='left', va='center')
    axs[1].axis('off')

    fig.savefig(Path(base_dir, 'f2-vaccine_profiles-full_7.png'))
# %%
