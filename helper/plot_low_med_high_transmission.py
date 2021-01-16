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
    
    if loc == 'RI':
        pop = 1059361.0
        pop_frac = np.array([0.105, 0.123, 0.140, 0.127, 0.124, 0.135, 0.120, 0.074, 0.052])
    elif loc == 'MA':
        pop = 6897212.0
        pop_frac = np.array([0.10586466, 0.12243686, 0.14498786, 0.13384234, 0.12230812, 0.14064248, 0.11801015, 0.06958116, 0.04232637])
    elif loc == 'PA':
        pop = 12800721.0
        pop_frac = np.array([0.11160395, 0.12229803, 0.13156525, 0.12581869, 0.11809624, 0.13878546, 0.12701660, 0.07657303, 0.04824275])
    
    savefig_dir = Path(cpp_dir.parent, 'sim_output/20210111-figures')
    out_dir_1 = Path(cpp_dir.parent, 'sim_output/20201226-low-allvac-{}-tot_300k'.format(loc))
    out_dir_2 = Path(cpp_dir.parent, 'sim_output/20201226-med-allvac-{}-tot_300k'.format(loc))
    out_dir_3 = Path(cpp_dir.parent, 'sim_output/20201226-high-allvac-{}-tot_300k'.format(loc))
    # out_dir_1 = Path(cpp_dir.parent, 'sim_output/20210109-low-allvac-{}-tot_1800k'.format(loc))
    # out_dir_2 = Path(cpp_dir.parent, 'sim_output/20210109-med-allvac-{}-tot_1800k'.format(loc))
    # out_dir_3 = Path(cpp_dir.parent, 'sim_output/20210109-high-allvac-{}-tot_1800k'.format(loc))    
    beta_1 = np.genfromtxt(Path(out_dir_1, '{}-low-beta0p1.txt'.format(loc)), delimiter=' ')
    beta_2 = np.genfromtxt(Path(out_dir_2, '{}-med-beta0p1.txt'.format(loc)), delimiter=' ')
    beta_3 = np.genfromtxt(Path(out_dir_3, '{}-high-beta0p1.txt'.format(loc)), delimiter=' ')
    
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
                  '16-29 only', '30-59 only', '60+ only',
                  '50/50 allocation among 20-39 and 60+', '25/75 allocation among 20-39 and 60+', '75/25 allocation among 20-39 and 60+', 
                  'random among 20-49 and 70+', '50/50 allocation among 20-49 and 70+', '25/75 allocation among 20-49 and 70+', '75/25 allocation among 20-49 and 70+']
    nsim = len(sim_labels)
    
    nat_imm_a=[540]
    normalcyday_a=[1098]
    vac_duration_a=[540]
    vac_halflife_a=[180, 180, 180, 360, 360, 540, 540]
    vac_slope_a =  [2,   3,   4,   1.5, 2,   1.5, 2]
    
    startday_1=370
    endday_1=429
    dpd_1=5000 if loc == 'RI' else 30000

    startday_2=370
    endday_2=429
    dpd_2=5000 if loc == 'RI' else 30000

    startday_3=370
    endday_3=429
    dpd_3=5000 if loc == 'RI' else 30000

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
    for out_dir, startday, endday, dpd in zip([out_dir_1, out_dir_2, out_dir_3], 
                                              [startday_1, startday_2, startday_3], 
                                              [endday_1, endday_2, endday_3], 
                                              [dpd_1, dpd_2, dpd_3]):
        if not Path(out_dir, 'AAA_s{}_e{}_d{}.dat'.format(startday, endday, dpd)).exists():
            ## first set of runs
            rundir = 's{}_e{}_d{}_i{}_c{}_u{}_t{}_p{}'.format(startday, endday, dpd, nat_imm_a[0], normalcyday_a[0], vac_duration_a[0], vac_halflife_a[0], vac_slope_a[0])
            # run_name = 'f426_t467_dpd10000' 
            A = Read_txt_To_ndarray(str(Path(out_dir, rundir)),
                                total_sim=nsim,
                                filename_format="run_output_{}.tdl", 
                                traj_ncol=int(sim.indices['DIMENSION']) + 1,
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
    with open(Path(out_dir_3, 'AAA_s{}_e{}_d{}.dat'.format(startday_3, endday_3, dpd_3)), 'rb') as f: AAA_3 = pcload(f)
    
    beta_1 = np.insert(beta_1, 0, [beta_1[0]]* (AAA_1.shape[-2] - beta_1.shape[0] +1))
    beta_2 = np.insert(beta_2, 0, [beta_2[0]]* (AAA_2.shape[-2] - beta_2.shape[0] +1))
    beta_3 = np.insert(beta_3, 0, [beta_3[0]]* (AAA_3.shape[-2] - beta_3.shape[0] +1))
    
    # %%
    ##      0. no vaccine
    # sim_labels = ['no vaccine', 'random', 
    #               '16-29 only', '30-59 only', '60+ only',
    #               '50/50 allocation among 20-39 and 60+', '25/75 allocation among 20-39 and 60+', '75/25 allocation among 20-39 and 60+', 
    #               'random among 20-49 and 70+', '50/50 allocation among 20-49 and 70+', '25/75 allocation among 20-49 and 70+', '75/25 allocation among 20-49 and 70+']and 70+', 'random 50/50 among 20-49 and 70+', 'random 25/75 among 20-49 and 70+', 'random 75/25 among 20-49 and 70+']
    # nat_imm_a=[540]
    # normalcyday_a=[1098]
    # vac_duration_a=[540]
    # vac_halflife_a=[180, 180, 180, 360, 360, 540, 540]
    # vac_slope_a =  [2,   3,   4,   1.5, 2,   1.5, 2]
    # %%
    subsim = [x for x in range(8)] + [x for x in range(9,nsim)] # [x for x in range(0,nsim)] # 
    subvac = [0] # [x for x in range(0, len(vac_halflife_a)) ]
    len_hl_a = len(subvac)
    norm = 0
    natimm = 0
    vacdur = 0
    ref = None ## the reference scenario, [0, nsim); if None, plot absolute values
    date_zero = pd.Timestamp('2019-12-31')
    date_startplt = pd.Timestamp('2020-08-15') # pd.Timestamp('2020-12-15') # 
    date_cutoff = pd.Timestamp('2021-06-30') # pd.Timestamp('2021-06-30') # 
    date_startvac = pd.to_datetime(min(startday_1, startday_2), unit='D', origin=date_zero)

    if date_startplt > date_zero and (date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        startrow_idx = int((date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) 
    else:
        startrow_idx = 0
    if date_cutoff > date_startplt and (date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        stoprow_idx = int((date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) + 1
    else:
        stoprow_idx = AAA_1.shape[-2]
        date_cutoff = pd.to_datetime(int(AAA_1[0,0,0,0, 0,-1,0]), unit='D', origin=date_zero )
    
    # %%
    attackrate_rightshift = 6 ## right-shift the J viewing window (i.e. shifting the whole curve left)
    seroprev_leftshift = 21 ## left-shift the J viewing window (i.e. shifting the whole curve right)

    ## seroprev
    asymp_frac_davies = np.array([0.71, 0.79, 0.73, 0.67, 0.60, 0.51, 0.37, 0.31, 0.31])
    asymp_frac_simple = np.array([0.381, 0.381, 0.460, 0.460, 0.431, 0.431, 0.409, 0.409, 0.130])
    asymp_frac = asymp_frac_davies # asymp_frac_simple


    ## sero-prevalence from date_startplt to date_cutoff
    dat_jz_ages = AAA_1[natimm, norm, vacdur, subvac[0], :, (startrow_idx-seroprev_leftshift):(stoprow_idx+attackrate_rightshift) , (sim.indices['STARTJ']+1):(sim.indices['STARTJ']+1+sim.indices['NUMAC']) ] 
    dat_jz_ages_w_asymp = dat_jz_ages / (1.0 - asymp_frac)
    
    attackrate_1 = np.sum(dat_jz_ages_w_asymp[:,:,:], axis=-1 ) / pop
    attackrate_1 = attackrate_1[:, (seroprev_leftshift+attackrate_rightshift):]

    # dat_jz_ages_w_asymp += AAA_1[natimm, norm, vacdur, subvac[0], :,  startrow_idx:stoprow_idx , (sim.indices['STARTJZ_1']+1):(sim.indices['STARTJZ_1']+1+sim.indices['NUMAC']) ]
    for nz in range(int(sim.indices['NUMZ_1'])):
        dat_jz_ages_w_asymp += AAA_1[natimm, norm, vacdur, subvac[0], :,  (startrow_idx-seroprev_leftshift):(stoprow_idx+attackrate_rightshift) , (sim.indices['STARTZ_1']+1+sim.indices['NUMAC']*nz):(sim.indices['STARTZ_1']+1+sim.indices['NUMAC']*(nz+1)) ]

    seroprev_1 = np.sum(dat_jz_ages_w_asymp[:,:,:], axis=-1 ) / pop
    seroprev_1 = seroprev_1[:, :-(seroprev_leftshift+attackrate_rightshift)]

    ###
    dat_jz_ages = AAA_2[natimm, norm, vacdur, subvac[0], :, (startrow_idx-seroprev_leftshift):(stoprow_idx+attackrate_rightshift) , (sim.indices['STARTJ']+1):(sim.indices['STARTJ']+1+sim.indices['NUMAC']) ] 
    dat_jz_ages_w_asymp = dat_jz_ages / (1.0 - asymp_frac)

    attackrate_2 = np.sum(dat_jz_ages_w_asymp[:,:,:], axis=-1 ) / pop
    attackrate_2 = attackrate_2[:, (seroprev_leftshift+attackrate_rightshift):]

    # dat_jz_ages_w_asymp += AAA_2[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , (sim.indices['STARTJZ_1']+1):(sim.indices['STARTJZ_1']+1+sim.indices['NUMAC']) ]
    for nz in range(int(sim.indices['NUMZ_1'])):
        dat_jz_ages_w_asymp += AAA_2[natimm, norm, vacdur, subvac[0], :,  (startrow_idx-seroprev_leftshift):(stoprow_idx+attackrate_rightshift) , (sim.indices['STARTZ_1']+1+sim.indices['NUMAC']*nz):(sim.indices['STARTZ_1']+1+sim.indices['NUMAC']*(nz+1)) ]

    seroprev_2 = np.sum(dat_jz_ages_w_asymp[:,:,:], axis=-1 ) / pop
    seroprev_2 = seroprev_2[:, :-(seroprev_leftshift+attackrate_rightshift)]

    ###
    dat_jz_ages = AAA_3[natimm, norm, vacdur, subvac[0], :, (startrow_idx-seroprev_leftshift):(stoprow_idx+attackrate_rightshift) , (sim.indices['STARTJ']+1):(sim.indices['STARTJ']+1+sim.indices['NUMAC']) ] 
    dat_jz_ages_w_asymp = dat_jz_ages / (1.0 - asymp_frac)

    attackrate_3 = np.sum(dat_jz_ages_w_asymp[:,:,:], axis=-1 ) / pop
    attackrate_3 = attackrate_3[:, (seroprev_leftshift+attackrate_rightshift):]

    # dat_jz_ages_w_asymp += AAA_3[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , (sim.indices['STARTJZ_1']+1):(sim.indices['STARTJZ_1']+1+sim.indices['NUMAC']) ]
    for nz in range(int(sim.indices['NUMZ_1'])):
        dat_jz_ages_w_asymp += AAA_3[natimm, norm, vacdur, subvac[0], :,  (startrow_idx-seroprev_leftshift):(stoprow_idx+attackrate_rightshift) , (sim.indices['STARTZ_1']+1+sim.indices['NUMAC']*nz):(sim.indices['STARTZ_1']+1+sim.indices['NUMAC']*(nz+1)) ]

    seroprev_3 = np.sum(dat_jz_ages_w_asymp[:,:,:], axis=-1 ) / pop
    seroprev_3 = seroprev_3[:, :-(seroprev_leftshift+attackrate_rightshift)]

    ## cumulative deaths, S class, I class
    D_1 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('D', traj=A) + sim.Get_Compartment_All_Ages('DHOSP', traj=A),
                                        AAA_1[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    D_2 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('D', traj=A) + sim.Get_Compartment_All_Ages('DHOSP', traj=A),
                                        AAA_2[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    D_3 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('D', traj=A) + sim.Get_Compartment_All_Ages('DHOSP', traj=A),
                                        AAA_3[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    I_1 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('I', traj=A),
                                        AAA_1[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    I_2 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('I', traj=A),
                                        AAA_2[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    I_3 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('I', traj=A),
                                        AAA_3[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] ))) 
    S_1 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('S', traj=A),
                                        AAA_1[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    S_2 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('S', traj=A),
                                        AAA_2[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))
    S_3 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('S', traj=A),
                                        AAA_3[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))

    J_1 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('J', traj=A),
                                        AAA_1[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    J_2 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('J', traj=A),
                                        AAA_2[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] )))  
    J_3 = np.array(list(map(lambda A: sim.Get_Compartment_All_Ages('J', traj=A),
                                        AAA_3[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx , : ] ))) 
    
    # %%
    # xdates = [pd.Timestamp('2020-09-01'), pd.Timestamp('2020-10-01'), pd.Timestamp('2020-11-01'), pd.Timestamp('2020-12-01'),
    #          pd.Timestamp('2021-01-01'), pd.Timestamp('2021-02-01'), pd.Timestamp('2021-03-01'), pd.Timestamp('2021-04-01'),
    #          pd.Timestamp('2021-05-01'), pd.Timestamp('2021-06-01')]
    xdates = [pd.Timestamp('2020-10-31'), pd.Timestamp('2020-11-30'), pd.Timestamp('2020-12-31'), 
              pd.Timestamp('2021-06-30')]
    for x in xdates:
        print(x.strftime('%Y-%m-%d'))
        print(attackrate_2[:, (x-date_startplt).days])
        print(seroprev_2[:, (x-date_startplt).days ])
        print(J_2[:, (x-date_startplt).days])

    # %%
    # cl_mix = ['mediumblue', 'xkcd:goldenrod', 'darkgreen', 'crimson', 'lime', 'orchid'] * 6
    # ls_mix = ['-'] * 6 + ['--'] * 6 +  [':'] * 6 + [(0,(1,5,1,5,3,5))] * 6 + [(0,(1,5,1,5,3,5,3,5))] * 6 + [(0, (1, 10))] * 6
    cl_mix = ['dimgrey', 'dimgrey', 
              'orchid', 'darkturquoise', 'seagreen', 
              'mediumblue', 'darkgreen', 'orangered', 
              'blue', 'green', 'red']
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
    # x_rightshift = -6
    fig, axs = plt.subplots(nrows=6, ncols=3, sharex='col', sharey='row', figsize=(27,24),
                            gridspec_kw={'left': 0.07, 'right': 0.70, 'bottom':0.12, 'top': 0.95,
                                         'wspace':0.04, 'hspace':0.11})
    x = AAA_1[0, 0, 0, 0, 0, startrow_idx:stoprow_idx , 0 ]
    for ax, srp, ttl in zip(axs[0,:].flat, [beta_1, beta_2, beta_3], ['low transmission', 'medium transmission', 'high transmission']):
        ax.plot(x, srp[startrow_idx:stoprow_idx], lw=4, color='black', ls='-')
        ax.set_title(ttl, fontsize=26)
        ax = Set_Axes_Xticks(ax, x[0], x[-1], tick_interval=30)

    for ax, srp in zip(axs[1,:].flat, [S_1, S_2, S_3]):
        for iss, s in enumerate(subsim):
            ax.plot(x, srp[s,:], label=sim_labels[s], lw=3, alpha=0.6)
        ax = Set_Axes_Xticks(ax, x[0], x[-1], tick_interval=30)

    for ax, srp in zip(axs[2,:].flat, [I_1, I_2, I_3]):
        for iss, s in enumerate(subsim):
            ax.plot(x, srp[s,:], label=sim_labels[s], lw=3, alpha=0.6)
        ax = Set_Axes_Xticks(ax, x[0], x[-1], tick_interval=30)
        if loc == 'RI':
            ax.set_yticks(range(0, 23001, 4000))
        if loc == 'MA':
            ax.set_yticks(range(0, 230001, 40000))

    for ax, srp in zip(axs[3,:].flat, [D_1, D_2, D_3]):
        for iss, s in enumerate(subsim):
            ax.plot(x, srp[s,:], label=sim_labels[s], lw=3, alpha=0.6)
        ax = Set_Axes_Xticks(ax, x[0], x[-1], tick_interval=30)
    
    for ax, srp in zip(axs[4,:].flat, [attackrate_1, attackrate_2, attackrate_3]):
        for iss, s in enumerate(subsim):
            ax.plot(x, srp[s,:], label=sim_labels[s], lw=3, alpha=0.6)
        ax = Set_Axes_Xticks(ax, x[0], x[-1], tick_interval=30)
    
    for ax, srp in zip(axs[5,:].flat, [seroprev_1, seroprev_2, seroprev_3]):
        for iss, s in enumerate(subsim):
            ax.plot(x, srp[s,:], label=sim_labels[s], lw=3, alpha=0.6)
        ax = Set_Axes_Xticks(ax, x[0], x[-1], tick_interval=30)
    
    axs[0,0].set_ylabel('Population mixing\n parameter', fontsize=25)
    axs[1,0].set_ylabel('S class', fontsize=25)
    axs[2,0].set_ylabel('I class', fontsize=25)
    axs[3,0].set_ylabel('Cumulative deaths', fontsize=25)
    axs[4,0].set_ylabel('Attack rate', fontsize=25)
    axs[5,0].set_ylabel('Seroprevalence', fontsize=25)
    axs[1,-1].legend(loc='upper left', bbox_to_anchor=(1.03, 1.0), fontsize=25)

    xdates = [pd.Timestamp('2020-09-01'), pd.Timestamp('2020-10-01'), pd.Timestamp('2020-11-01'), pd.Timestamp('2020-12-01'),
             pd.Timestamp('2021-01-01'), pd.Timestamp('2021-02-01'), pd.Timestamp('2021-03-01'), pd.Timestamp('2021-04-01'),
             pd.Timestamp('2021-05-01'), pd.Timestamp('2021-06-01'), pd.Timestamp('2021-07-01')]
    xdates_excl = set([pd.Timestamp('2020-10-01'), pd.Timestamp('2020-12-01'),
                   pd.Timestamp('2021-02-01'), pd.Timestamp('2021-04-01'), pd.Timestamp('2021-06-01')])
    xdates_ticks = [ (x-date_zero).days for x in xdates ]
    xdates_labels = [ x.strftime('%b %d') for x in xdates ]
    for ixt, xt in enumerate(xdates):
        if xt in xdates_excl:
            xdates_labels[ixt] = ''
        if xt == pd.Timestamp('2021-01-01'):
            xdates_labels[ixt] = xt.strftime('%Y - %b %d')
    for ax in axs.flat:
        ax.set_xticks(xdates_ticks)
        ax.set_xticklabels(xdates_labels, fontsize=24)
        ax.axvline(x=(pd.Timestamp('2021-01-01') - date_zero).days, color='black', ls='--', lw=2)
        ax.margins(x=0.07, y=0.1)


    fig.suptitle('{}, vaccine campaign starts on {}, ends on {} or {}, with {} doses per day'.format(loc, 
                date_startvac.strftime('%d %b %Y'), 
                pd.to_datetime(endday_1, unit='D', origin=date_zero).strftime('%d %b %Y'),
                pd.to_datetime(endday_2, unit='D', origin=date_zero).strftime('%d %b %Y'),
                dpd_1 ), y=0.02)

    fig.savefig(Path(savefig_dir, '{}-low_med_high_transmission-i{}_c{}_u{}_e{}_d{}.png'.format(loc, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur], endday_1, dpd_1)))





    
    # %%
    ## outcomes summary
    for AAA, d_sp, atr, serop, fn in zip([AAA_1, AAA_2, AAA_3], [endday_1, endday_2, endday_3], 
                                    [attackrate_1, attackrate_2, attackrate_3],
                                    [seroprev_1, seroprev_2, seroprev_3],
                                    ['low-tot_1800k', 'med-tot_1800k', 'high-tot_1800k']):
        ## sero-prevalence
        date_endvac = pd.to_datetime(d_sp, unit='D', origin=date_zero)
        # row_idx = (date_endvac - pd.Timestamp('2020-03-01')).days
        # row_idx = (date_endvac - date_zero).days - int(AAA[0, 0, 0, 0, 0, 0 , 0] )

        # dat_jz_ages = AAA[natimm, norm, vacdur, subvac[0], :, (row_idx-1):(row_idx+1) , (sim.indices['STARTJ']+1):(sim.indices['STARTJ']+1+sim.indices['NUMAC']) ] 
        # dat_jz_ages_w_asymp = dat_jz_ages / (1.0 - asymp_frac)
        # dat_jz_ages_w_asymp += AAA[natimm, norm, vacdur, subvac[0], :, (row_idx-1):(row_idx+1) , (sim.indices['STARTJZ_1']+1):(sim.indices['STARTJZ_1']+1+sim.indices['NUMAC']) ]

        # serop = np.sum(dat_jz_ages_w_asymp[:,-1,:], axis=-1 ) / pop

        f = pd.DataFrame()
        f['strategy'] = sim_labels
        f['cumulative_cases_by_{}'.format(date_cutoff.strftime("%Y%m%d"))] = list(map(lambda A: sim.Get_Compartment_All_Ages("J", traj=A)[-1], 
                                                                                        AAA[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx,:] ))
        f['cumulative_hospitalizations_by_{}'.format(date_cutoff.strftime("%Y%m%d"))] = list(map(lambda A: sim.Get_Compartment_All_Ages("K", traj=A)[-1], 
                                                                                        AAA[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx,:] ))
        f['cumulative_deaths_by_{}'.format(date_cutoff.strftime("%Y%m%d"))] = list(map(lambda A: (sim.Get_Compartment_All_Ages("D", traj=A) + sim.Get_Compartment_All_Ages("DHOSP", traj=A))[-1], 
                                                                                        AAA[natimm, norm, vacdur, subvac[0], :, startrow_idx:stoprow_idx,:] ))
        # f['attackrate_by_{}'.format(date_endvac.strftime("%Y%m%d"))] = atr[:,-1]
        # f['seroprev_by_{}'.format(date_endvac.strftime("%Y%m%d"))] = serop[:,-1]
        f['attackrate_by_{}'.format(date_cutoff.strftime("%Y%m%d"))] = atr[:,-1]
        f['seroprev_by_{}'.format(date_cutoff.strftime("%Y%m%d"))] = serop[:,-1]

        f.to_csv(Path(savefig_dir, '{}-{}-summary_outcomes_{}.csv'.format(loc, fn, date_cutoff.strftime("%Y%m%d")) ), index=False)

    
# %%
