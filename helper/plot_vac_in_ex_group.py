#!/usr/bin/env python3
# %%
import covid_sim
from plot_helper import Set_Axes_Xticks
from gather_output import Read_txt_To_ndarray
from pathlib import Path
# from pickle import load as pcload
from _pickle import load as pcload ## cPickle
from _pickle import dump as pcdump
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
# from pymc3.stats import hpd
# from matplotlib.gridspec import GridSpec ## for figure layour
from matplotlib import rc_file

# plt.rcdefaults()
rc_file(Path(Path.home(), 'Code/covid19-vaccine/helper/matplotlibrc-custom'))


# %%
if __name__ == "__main__":
    # %%
    cpp_dir = Path(Path.home(), 'Code/covid19-vaccine/cpp-v6-test-vaccination/') 
    sim = covid_sim.COVID_SIM(cpp_dir)
    numac = int(sim.indices['NUMAC'])

    loc = 'MA'
    savefig_dir = Path(cpp_dir.parent, 'sim_output/20210111-figures')
    out_dir = Path(cpp_dir.parent, 'sim_output/20210111-med-inexgrp-{}-tot_1800k/'.format(loc))
    conf = pd.read_csv(Path(out_dir, 'config-in_ex_group.csv'))
    nsim = conf.shape[0]
    startday = 370
    endday = 429
    dpd = 5000 if loc == 'RI' else 30000
    conf_sim = '{}, startday {}, endday {}, dpd {}, nat_immu 540, vac_protect 540, vac_halflife 360, vac_slope 2'.format(loc, startday, endday, dpd)
    
    # %%
    if not Path(out_dir, 'AAA-in_ex_group-s{}_e{}_d{}.dat'.format(startday, endday, dpd)).exists():
        ## read all sim output
        AAA = Read_txt_To_ndarray(str(Path(out_dir, 's{}_e{}_d{}'.format(startday, endday, dpd))),
                            total_sim=nsim,
                            filename_format="run_output_{}.tdl", 
                            traj_ncol=int(sim.indices['DIMENSION']) +1,
                            binary_output=True,
                            death_rate_only=False)

        ## save msater array for later loading
        with open(Path(out_dir, 'AAA-in_ex_group-s{}_e{}_d{}.dat'.format(startday, endday, dpd)), 'wb') as f: pcdump(AAA, f) 
    # %%
    with open(Path(out_dir, 'AAA-in_ex_group-s{}_e{}_d{}.dat'.format(startday, endday, dpd)), 'rb') as f: AAA = pcload(f)
    
    # %%
    ## prepare sim indices
    # incl_00_idx = conf.loc[ conf['incl_00'] == 1 , 'sim'].to_list()
    # incl_10_idx = conf.loc[ conf['incl_10'] == 1 , 'sim'].to_list()
    # incl_20_idx = conf.loc[ conf['incl_20'] == 1 , 'sim'].to_list()
    # incl_30_idx = conf.loc[ conf['incl_30'] == 1 , 'sim'].to_list()
    # incl_40_idx = conf.loc[ conf['incl_40'] == 1 , 'sim'].to_list()
    # incl_50_idx = conf.loc[ conf['incl_50'] == 1 , 'sim'].to_list()
    # incl_60_idx = conf.loc[ conf['incl_60'] == 1 , 'sim'].to_list()
    # incl_70_idx = conf.loc[ conf['incl_70'] == 1 , 'sim'].to_list()
    # incl_80_idx = conf.loc[ conf['incl_80'] == 1 , 'sim'].to_list()
    # excl_00_idx = conf.loc[ conf['incl_00'] == 0 , 'sim'].to_list()
    # excl_10_idx = conf.loc[ conf['incl_10'] == 0 , 'sim'].to_list()
    # excl_20_idx = conf.loc[ conf['incl_20'] == 0 , 'sim'].to_list()
    # excl_30_idx = conf.loc[ conf['incl_30'] == 0 , 'sim'].to_list()
    # excl_40_idx = conf.loc[ conf['incl_40'] == 0 , 'sim'].to_list()
    # excl_50_idx = conf.loc[ conf['incl_50'] == 0 , 'sim'].to_list()
    # excl_60_idx = conf.loc[ conf['incl_60'] == 0 , 'sim'].to_list()
    # excl_70_idx = conf.loc[ conf['incl_70'] == 0 , 'sim'].to_list()
    # excl_80_idx = conf.loc[ conf['incl_80'] == 0 , 'sim'].to_list()

    # %%
    ## summarizing cases, hospitalizations, deaths for all ages
    AAA_J = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages("J",traj=A), AAA ) ))
    AAA_K = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages("K",traj=A), AAA ) ))
    AAA_D = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages("D",traj=A) + sim.Get_Compartment_All_Ages("DHOSP",traj=A), AAA ) ))
    # AAA_JZ_1 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages("JZ_1",traj=A), AAA ) ))
    # np.median(AAA_JZ_1[1:,-1])

    # %%
    ## prepare array for violin plot
    date_zero = pd.Timestamp('2019-12-31')
    date_cutoff = pd.Timestamp('2021-06-30')
    if (date_cutoff - date_zero).days - AAA[0,0,0] < AAA.shape[1]:
        row_idx = int((date_cutoff - date_zero).days - AAA[0,0,0])
    else:
        row_idx = -1
    ##
    ## plt_J, plt_K, plt_D
    ## shape: (nsim/2, NUMAC*2)
    ## the first NUMAC columns are for `inclusion` of a certain age group (each column is one age group)
    ## the last NUMAC columns are for `exclusion` of a certain age group
    plt_J = np.empty(( nsim//2, numac*2))
    plt_K = np.empty(( nsim//2, numac*2))
    plt_D = np.empty(( nsim//2, numac*2))
    for i in range(numac):
        if i == 0: 
            incl_idx = conf.loc[ conf['incl_00'] == 1 , 'sim'].to_list()
            excl_idx = conf.loc[ conf['incl_00'] == 0 , 'sim'].to_list()
        else: 
            incl_idx = conf.loc[ conf['incl_{}'.format(i*10)] == 1 , 'sim'].to_list()
            excl_idx = conf.loc[ conf['incl_{}'.format(i*10)] == 0 , 'sim'].to_list()
        plt_J[:,i] = AAA_J[incl_idx, row_idx]
        plt_K[:,i] = AAA_K[incl_idx, row_idx]
        plt_D[:,i] = AAA_D[incl_idx, row_idx]
        plt_J[:,numac+i] = AAA_J[excl_idx, row_idx]
        plt_K[:,numac+i] = AAA_K[excl_idx, row_idx]
        plt_D[:,numac+i] = AAA_D[excl_idx, row_idx]
    

    # %%
    ## export min, max, median, IQR to csv
    plt_stats = np.array([list(np.min(plt_J[: ,:numac], axis=0)) ] + \
                         [list(np.max(plt_J[: ,:numac], axis=0)) ] + \
                         [x for x in np.quantile( plt_J[: ,:numac] , [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_J[1:,numac:], axis=0)) ] + \
                         [list(np.max(plt_J[1:,numac:], axis=0)) ] + \
                         [x for x in np.quantile( plt_J[1:,numac:], [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_K[: ,:numac], axis=0)) ] + \
                         [list(np.max(plt_K[: ,:numac], axis=0)) ] + \
                         [x for x in np.quantile( plt_K[: ,:numac] , [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_K[1:,numac:], axis=0)) ] + \
                         [list(np.max(plt_K[1:,numac:], axis=0)) ] + \
                         [x for x in np.quantile( plt_K[1:,numac:], [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_D[: ,:numac], axis=0)) ] + \
                         [list(np.max(plt_D[: ,:numac], axis=0)) ] + \
                         [x for x in np.quantile( plt_D[: ,:numac] , [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_D[1:,numac:], axis=0)) ] + \
                         [list(np.max(plt_D[1:,numac:], axis=0)) ] + \
                         [x for x in np.quantile( plt_D[1:,numac:], [0.5, 0.25, 0.75], axis=0)])

    f = pd.DataFrame({'cumulative': np.repeat(['cases', 'hospitalizations', 'deaths'], 10),
                      'include': np.tile( np.repeat(['yes', 'no'], 5), 3 ),
                      'stats': np.tile(['min', 'max', 'q0.5', 'q0.25', 'q0.75'], 6),
                      '0_9': plt_stats[:,0],
                      '10_19': plt_stats[:,1], 
                      '20_29': plt_stats[:,2], 
                      '30_39': plt_stats[:,3], 
                      '40_49': plt_stats[:,4], 
                      '50_59': plt_stats[:,5],
                      '60_69': plt_stats[:,6], 
                      '70_79': plt_stats[:,7],
                      '80+': plt_stats[:,8] })

    # f.to_csv(Path(out_dir, '{}-med_beta-incl_excl-s{}_e{}_d{}-{}-summary-JKD.csv'.format(loc, startday, endday, dpd, date_cutoff.strftime("%y%m%d"))), index=False)
    f.to_csv(Path(savefig_dir, '{}-med_beta-incl_excl-summary-JKD.csv'.format(loc)), index=False)

    # %%
    ## export %change (if included) in min, max, median, IQR to csv
    plt_J_delta = (plt_J[:,:numac] - plt_J[:,numac:])/ plt_J[:,numac:]
    plt_K_delta = (plt_K[:,:numac] - plt_K[:,numac:])/ plt_K[:,numac:]
    plt_D_delta = (plt_D[:,:numac] - plt_D[:,numac:])/ plt_D[:,numac:]

    plt_stats = np.array([list(np.min(plt_J_delta[1:,:], axis=0)) ] + \
                         [list(np.max(plt_J_delta[1:,:], axis=0)) ] + \
                         [x for x in np.quantile( plt_J_delta[1:,:] , [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_K_delta[1:,:], axis=0)) ] + \
                         [list(np.max(plt_K_delta[1:,:], axis=0)) ] + \
                         [x for x in np.quantile( plt_K_delta[1:,:] , [0.5, 0.25, 0.75], axis=0)] + \
                         [list(np.min(plt_D_delta[1:,:], axis=0)) ] + \
                         [list(np.max(plt_D_delta[1:,:], axis=0)) ] + \
                         [x for x in np.quantile( plt_D_delta[1:,:] , [0.5, 0.25, 0.75], axis=0)] )

    f = pd.DataFrame({'cumulative': np.repeat(['cases', 'hospitalizations', 'deaths'], 5),
                      'include_over_exclude': np.tile( np.repeat(['yes'], 5), 3 ),
                      'stats': np.tile(['min', 'max', 'q0.5', 'q0.25', 'q0.75'], 3),
                      '0_9': plt_stats[:,0],
                      '10_19': plt_stats[:,1], 
                      '20_29': plt_stats[:,2], 
                      '30_39': plt_stats[:,3], 
                      '40_49': plt_stats[:,4], 
                      '50_59': plt_stats[:,5],
                      '60_69': plt_stats[:,6], 
                      '70_79': plt_stats[:,7],
                      '80+': plt_stats[:,8] })

    # f.to_csv(Path(out_dir, '{}-med_beta-incl_excl-s{}_e{}_d{}-{}-summary_incl_over_excl-JKD.csv'.format(loc, startday, endday, dpd, date_cutoff.strftime("%y%m%d"))), index=False)
    f.to_csv(Path(savefig_dir, '{}-med_beta-incl_excl-{}-summary_incl_over_excl-JKD.csv'.format(loc, date_cutoff.strftime("%y%m%d"))), index=False)
    
    # %%
    ## violin plot
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(numac*2, 18), 
                            gridspec_kw={'top':0.95, 'bottom':0.07, 'hspace':0.4})
    for ax, plt_a, ttl in zip(axs.flat, [plt_J, plt_K, plt_D], ['cases', 'hospitalizations', 'deaths']):
        ax.violinplot(plt_a[:, :numac], positions=[i-0.3 for i in np.linspace(0,14,numac)])
        ax.violinplot(plt_a[1:, numac:], positions=[i+0.3 for i in np.linspace(0,14,numac)])
        ax.set_xticks(np.linspace(0,14,numac))
        ax.set_xticklabels(['0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'])
        if row_idx == -1:
            ax.set_title('Cumulative {} by the end of the simulation'.format( ttl))
        else:
            ax.set_title('Cumulative {} by {}'.format(ttl, date_cutoff.strftime('%d %B %Y')))
    fig.suptitle(conf_sim, y=0.02)

    # fig.savefig(Path(out_dir, '{}-med_beta-incl_excl-s{}_e{}_d{}-{}_no_baseline.png'.format(loc, startday, endday, dpd, date_cutoff.strftime("%y%m%d"))))
    fig.savefig(Path(savefig_dir, '{}-med_beta-incl_excl-no_baseline.png'.format(loc)))
# %%
