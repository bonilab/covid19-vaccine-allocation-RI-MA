#!/usr/bin/env python3
# %%
import covid_sim
from plot_helper import Set_Axes_Xticks
from plot_efficacy import Get_Efficacy
from gather_output import Read_txt_To_ndarray
from pathlib import Path
# from pickle import load as pcload
from _pickle import load as pcload ## cPickle
from _pickle import dump as pcdump
import pandas as pd 
import numpy as np 
from cycler import cycler
import matplotlib.pyplot as plt 
from matplotlib import cm
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
    # numac = int(sim.indices['NUMAC'])
    
    loc = 'RI'## MODIFY HERE
    trset = 'high'
    
    savefig_dir = Path(cpp_dir.parent, 'sim_output/20201226-figures')
    out_dir_1 = Path(cpp_dir.parent, 'sim_output/20201226-{}-allvac-{}-tot_50k'.format(trset, loc))
    out_dir_2 = Path(cpp_dir.parent, 'sim_output/20201226-{}-allvac-{}-tot_300k'.format(trset, loc))
    
    if loc == 'RI':
        pop = 1059361.0
        pop_frac = np.array([0.105, 0.123, 0.140, 0.127, 0.124, 0.135, 0.120, 0.074, 0.052])
    elif loc == 'MA':
        pop = 6897212.0
        pop_frac = np.array([0.10586466, 0.12243686, 0.14498786, 0.13384234, 0.12230812, 0.14064248, 0.11801015, 0.06958116, 0.04232637])
    elif loc == 'PA':
        pop = 12800721.0
        pop_frac = np.array([0.11160395, 0.12229803, 0.13156525, 0.12581869, 0.11809624, 0.13878546, 0.12701660, 0.07657303, 0.04824275])
    
    ## vaccination stratgies:
    ##      0. no vaccine
    ##      1. random
    ##      2. 16-29 only
    ##      3. 30-59 only
    ##      4. 60+ only
    ##      5. random 50/50 among 20-39 and 60+
    ##      6. random 25/75 among 20-39 and 60+
    ##      7. random 75/25 among 20-39 and 60+
    ##      8. random (60/40) among 20-49 and 70+
    ##      9. random 50/50 among 20-49 and 70+
    ##     10. random 25/75 among 20-49 and 70+
    ##     11. random 75/25 among 20-49 and 70+
    sim_labels = ['no vaccine', 'random', 
                  '16-29', '30-59', '60+',
                  '50/50 allocation among 20-39 and 60+', '25/75 allocation among 20-39 and 60+', '75/25 allocation among 20-39 and 60+', 
                  'random among 20-49 and 70+', '50/50 allocation among 20-49 and 70+', '25/75 allocation among 20-49 and 70+', '75/25 allocation among 20-49 and 70+']
    nsim = len(sim_labels)
    
    nat_imm_a=[540]
    normalcyday_a=[1098]
    vac_duration_a=[540]
    vac_halflife_a=[180, 180, 180, 360, 360, 540, 540]
    vac_slope_a =  [2,   3,   4,   1.5, 2,   1.5, 2]
    
    startday_1=370
    endday_1=379
    dpd_1=5000 if loc == 'RI' else 30000

    startday_2=370
    endday_2=429
    dpd_2=5000 if loc == 'RI' else 30000

    ##
    ## AAA master array:
    #       natural_immunity
    #       normalcy_day 
    #       vaccine_protection_duration
    #       vaccine_halflife (paired with vaccine_slope)
    #       nsim 
    #       odesim_rows 
    #       odesim_columns
    # %%
    for out_dir, startday, endday, dpd in zip([out_dir_1, out_dir_2], [startday_1, startday_2], [endday_1, endday_2], [dpd_1, dpd_2]):
        if not Path(out_dir, 'AAA_s{}_e{}_d{}.dat'.format(startday, endday, dpd)).exists():
            ## first set of runs
            rundir = 's{}_e{}_d{}_i{}_c{}_u{}_t{}_p{}'.format(startday, endday, dpd, nat_imm_a[0], normalcyday_a[0], vac_duration_a[0], vac_halflife_a[0], vac_slope_a[0])
            # run_name = 'f426_t467_dpd10000' 
            A = Read_txt_To_ndarray(str(Path(out_dir, rundir)),
                                total_sim=nsim,
                                filename_format="run_output_{}.tdl", 
                                traj_ncol=int(sim.indices['DIMENSION']) +1,
                                binary_output=True,
                                death_rate_only=False)

            AAA = np.empty((len(nat_imm_a), len(normalcyday_a), len(vac_duration_a), len(vac_halflife_a), A.shape[0], A.shape[1], A.shape[2] ))
            for nii, ni in enumerate(nat_imm_a):
                for ndi, nd in enumerate(normalcyday_a):
                    for vdi, vd in enumerate(vac_duration_a):
                        for vhi, vh in enumerate(vac_halflife_a):
                            rundir = 's{}_e{}_d{}_i{}_c{}_u{}_t{}_p{}'.format(startday, endday, dpd, ni, nd, vd, vh, vac_slope_a[vhi])
                            AAA[nii, ndi, vdi, vhi, :,:,:] = Read_txt_To_ndarray(str(Path(out_dir, rundir)),
                                                                total_sim=nsim,
                                                                filename_format="run_output_{}.tdl", 
                                                                traj_ncol=int(sim.indices['DIMENSION']) +1,
                                                                binary_output=True,
                                                                death_rate_only=False)
            
            ## save msater array for later loading
            with open(Path(out_dir, 'AAA_s{}_e{}_d{}.dat'.format(startday, endday, dpd)), 'wb') as f: pcdump(AAA, f) 
    # %%
    with open(Path(out_dir_1, 'AAA_s{}_e{}_d{}.dat'.format(startday_1, endday_1, dpd_1)), 'rb') as f: AAA_1 = pcload(f)
    with open(Path(out_dir_2, 'AAA_s{}_e{}_d{}.dat'.format(startday_2, endday_2, dpd_2)), 'rb') as f: AAA_2 = pcload(f)
    
    
    # %%
    ##      0. no vaccine
    # sim_labels = ['no vaccine', 'random', 
    #               '16-29', '30-59', '60+',
    #               '50/50 allocation among 20-39 and 60+', '25/75 allocation among 20-39 and 60+', '75/25 allocation among 20-39 and 60+', 
    #               'random among 20-49 and 70+', '50/50 allocation among 20-49 and 70+', '25/75 allocation among 20-49 and 70+', '75/25 allocation among 20-49 and 70+']
    # nat_imm_a=[540]
    # normalcyday_a=[1098]
    # vac_duration_a=[540]
    # vac_halflife_a=[180, 180, 180, 360, 360, 540, 540]
    # vac_slope_a =  [2,   3,   4,   1.5, 2,   1.5, 2]
    # %%
    subsim = [x for x in range(8)] + [x for x in range(9,nsim)] # [x for x in range(0,nsim)] # 
    subvac = [5, 6] # [3, 4] # [0, 1, 2] #  [x for x in range(len(vac_halflife_a)) ]
    len_hl_a = len(subvac)
    len_subsim = len(subsim)
    norm = 0
    natimm = 0
    vacdur = 0
    ref = 1 ## the reference scenario, [0, nsim); if None, plot absolute values
    date_zero = pd.Timestamp('2019-12-31')
    date_startplt = pd.Timestamp('2020-12-25') # pd.Timestamp('2020-12-15') # 
    date_cutoff = pd.Timestamp('2021-06-30') # pd.Timestamp('2021-06-30') # 
    date_startvac = pd.to_datetime(min(startday_1, startday_2), unit='D', origin=date_zero)

    if date_startplt > date_zero and (date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        startrow_idx = int((date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) + 1
    else:
        startrow_idx = 0
    if date_cutoff > date_startplt and (date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        stoprow_idx = int((date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) + 2
    else:
        stoprow_idx = AAA_1.shape[-2]
    
    # %%
    # cl_mix = ['mediumblue', 'xkcd:goldenrod', 'darkgreen', 'crimson', 'lime', 'orchid'] * 6
    # ls_mix = ['-'] * 6 + ['--'] * 6 +  [':'] * 6 + [(0,(1,5,1,5,3,5))] * 6 + [(0,(1,5,1,5,3,5,3,5))] * 6 + [(0, (1, 10))] * 6
    cl_mix = ['dimgrey', 'dimgrey', 
                  'orchid', 'darkturquoise', 'seagreen', 
                  'mediumblue', 'green', 'orangered', 
                  'xkcd:blue', 'darkgreen', 'tab:red']
                #   'xkcd:goldenrod', 'blue', 'green', 'red', 
                #   'xkcd:goldenrod']
    ls_mix = [':','-',
              '-', '-', '-',
              '-.', '-.', '-.', 
              '--', '--', '--']
            #   '--', '--', '--', '--',
            #   '-']
    style_cycler = cycler(color=cl_mix, linestyle=ls_mix)
    # style_cycler = cycler(color=cl_mix[:nsim])

    plt.rc('axes', prop_cycle=style_cycler )
    
    # %%
    ## combine 10-day and 60-day campaign
    ##
    #### 3 subvac
    if len_hl_a == 3:
        fig, axs = plt.subplots(nrows=4, ncols=len_hl_a*2 +1, 
                                sharex='col', sharey='row',
                                figsize=(7*(len_hl_a*2 ), 6*4), dpi=120,
                                gridspec_kw={'hspace': 0.07, 'wspace': .16, 
                                            'bottom': 0.10, 'top': 0.93, 'left': 0.08, 'right': 0.81,
                                            'width_ratios': [20]*len_hl_a + [1] + [20]*len_hl_a})
    #### 2 subvac
    if len_hl_a == 2:
        fig, axs = plt.subplots(nrows=4, ncols=len_hl_a*2 +1, 
                                sharex='col', sharey='row',
                                figsize=(7.5*(len_hl_a*2 ), 6*4), dpi=120,
                                gridspec_kw={'hspace': 0.07, 'wspace': .16, 
                                            'bottom': 0.10, 'top': 0.93, 'left': 0.08, 'right': 0.73,
                                            'width_ratios': [20]*len_hl_a + [1] + [20]*len_hl_a})
        
    for ic, iv in enumerate(subvac + [-1] + subvac):
        if ic != len_hl_a:
            v = vac_halflife_a[iv]
            slope = vac_slope_a[iv]
            
            ## first row: efficacy curve
            eff_days = np.arange( (date_cutoff-date_startvac).days )
            eff = [Get_Efficacy(x, v, slope ) for x in eff_days ]
            axs[0,ic].plot((date_startvac - date_zero).days + eff_days , eff, color='black', ls='-', lw=4.5)
            axs[0,ic].set_title('vaccine halflife: {}\nvaccine slope: {}'.format(v, slope))
            
            AAA = AAA_1 if ic < len_hl_a else AAA_2
            ## second, third, fourth row: (relative) cumulative cases, hospitalization, deaths
            if ref is not None:
                base_J = sim.Get_Compartment_All_Ages("J",traj=AAA[natimm, norm, vacdur, iv, ref, startrow_idx:stoprow_idx,:])
                base_K = sim.Get_Compartment_All_Ages("K",traj=AAA[natimm, norm, vacdur, iv, ref, startrow_idx:stoprow_idx,:])
                base_D = sim.Get_Compartment_All_Ages("D",traj=AAA[natimm, norm, vacdur, iv, ref, startrow_idx:stoprow_idx,:]) + \
                        sim.Get_Compartment_All_Ages("DHOSP",traj=AAA[natimm, norm, vacdur,iv, ref, startrow_idx:stoprow_idx,:])
            else:
                base_J = base_K = base_D = 1
            for iss, s in enumerate(subsim):
                A = AAA[natimm, norm, vacdur, iv, s,startrow_idx:stoprow_idx,:]
                axs[1,ic].plot(A[:,0], (sim.Get_Compartment_All_Ages("J",traj=A)) / base_J, label=sim_labels[s], alpha=0.6, lw=4.5)
                axs[2,ic].plot(A[:,0], (sim.Get_Compartment_All_Ages("K",traj=A)) / base_K, label=sim_labels[s], alpha=0.6, lw=4.5)
                axs[3,ic].plot(A[:,0], (sim.Get_Compartment_All_Ages("D",traj=A) + sim.Get_Compartment_All_Ages("DHOSP",traj=A)) / base_D, label=sim_labels[s], alpha=0.6, lw=4.5)
                ## plot cumulative vaccinations
                # axs[4,ic].plot(A[:,0], (sim.Get_Compartment_All_Ages("JZ_1",traj=A)), label=sim_labels[s], alpha=0.6, lw=4.0)
    
    xdates = [pd.Timestamp('2021-01-01'), pd.Timestamp('2021-02-01'), pd.Timestamp('2021-03-01'), pd.Timestamp('2021-04-01'),
             pd.Timestamp('2021-05-01'), pd.Timestamp('2021-06-01'), pd.Timestamp('2021-07-01')]
    xdates_excl = set([pd.Timestamp('2021-02-01'), pd.Timestamp('2021-04-01'), pd.Timestamp('2021-06-01')])
    xdates_ticks = [ (x-date_zero).days for x in xdates ]
    xdates_labels = [ x.strftime('%b %d') for x in xdates ]
    for ixt, xt in enumerate(xdates):
        if xt in xdates_excl:
            xdates_labels[ixt] = ''
        if xt == pd.Timestamp('2021-01-01'):
            xdates_labels[ixt] = xt.strftime('%Y - %b %d')
    for ax in axs.flat:
        Set_Axes_Xticks(ax, (pd.Timestamp('2021-01-01') - date_zero).days, (date_cutoff - date_zero).days, tick_interval=30)
        ax.set_xticks(xdates_ticks)
        ax.set_xticklabels(xdates_labels, fontsize=24)
        ax.set_xlim(left=(date_startplt - date_zero).days, right=(date_cutoff - date_zero).days+1 )
        ax.axvline(x=(date_startvac- date_zero).days, color='grey', ls='--', lw=2)
    
    for ax in axs[0,:].flat:
        ax.set_yticks(np.arange(0,1.01,0.2))
        ax.set_yticklabels(['{:.0%}'.format(x) for x in  np.arange(0, 1.01, 0.2)])
    axs[0,0].set_ylabel('Vaccince efficacy')

    axs[1,-1].legend(loc='upper left', bbox_to_anchor=(1.03, 1.0), fontsize=24)

    if ref is not None:
        axs[1,0].set_ylabel('Cumulative cases\nrelative to the\n"{}" strategy'.format(sim_labels[ref]))
        axs[2,0].set_ylabel('Cumulative hospitalizations\nrelative to the\n"{}" strategy'.format(sim_labels[ref]))
        axs[3,0].set_ylabel('Cumulative deaths\nrelative to the\n"{}" strategy'.format(sim_labels[ref]))
    else:
        axs[1,0].set_ylabel('Cumulative cases')
        axs[2,0].set_ylabel('Cumulative hospitalizations')
        axs[3,0].set_ylabel('Cumulative deaths')

    for ax in axs[:, len_hl_a]:
        ax.axis('off')
    
    fig.suptitle('{}, vaccine campaign starts on {}, ends on {} (low supply) or on {} (ample supply), with {} doses per day'.format(loc, 
                        date_startvac.strftime('%d %b %Y'), 
                        pd.to_datetime(endday_1, unit='D', origin=date_zero).strftime('%d %b %Y'),
                        pd.to_datetime(endday_2, unit='D', origin=date_zero).strftime('%d %b %Y'),
                        dpd_1 ), y=0.02)

    
    if ref is not None: fig.savefig(Path(savefig_dir, '{}-{}_beta-allvac-t{}-combine-ref{}.png'.format(loc, trset, vac_halflife_a[subvac[0]] , ref )))
    else:               fig.savefig(Path(savefig_dir, '{}-{}_beta-allvac-t{}-combine-abs.png'.format(loc, trset, vac_halflife_a[subvac[0]])))



# %%
