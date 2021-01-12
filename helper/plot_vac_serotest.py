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

    loc = "RI"
    trset = 'med'
    savefig_dir = Path(cpp_dir.parent, 'sim_output/20210111-figures')
    out_dir_1 = Path(cpp_dir.parent, 'sim_output/20210111-{}-allvac-{}-tot_50k'.format(trset, loc))
    out_dir_2 = Path(cpp_dir.parent, 'sim_output/20210111-{}-allvac-{}-tot_50k-notest'.format(trset, loc))
    out_dir_3 = Path(cpp_dir.parent, 'sim_output/20210111-{}-allvac-{}-tot_300k'.format(trset, loc))
    out_dir_4 = Path(cpp_dir.parent, 'sim_output/20210111-{}-allvac-{}-tot_300k-notest'.format(trset, loc))
    
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
    endday_2=379
    dpd_2=5000 if loc == 'RI' else 30000

    startday_3=370
    endday_3=429 
    dpd_3=5000 if loc == 'RI' else 30000

    startday_4=370
    endday_4=429 
    dpd_4=5000 if loc == 'RI' else 30000

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
    for out_dir, startday, endday, dpd in zip([out_dir_1, out_dir_2, out_dir_3, out_dir_4], [startday_1, startday_2, startday_3, startday_4], [endday_1, endday_2, endday_3, endday_4], [dpd_1, dpd_2, dpd_3, dpd_4]):
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
    with open(Path(out_dir_3, 'AAA_s{}_e{}_d{}.dat'.format(startday_3, endday_3, dpd_3)), 'rb') as f: AAA_3 = pcload(f)
    with open(Path(out_dir_4, 'AAA_s{}_e{}_d{}.dat'.format(startday_4, endday_4, dpd_4)), 'rb') as f: AAA_4 = pcload(f)
        
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
    subsim = [10, 11] # [x for x in range(1,nsim)] # [1] + [x for x in range(3,nsim)] #
    subsim = sorted(subsim)
    subvac = [0, 1, 2] # [x for x in range(0, len(vac_halflife_a)) ]
    subvac = sorted(subvac)
    norm = 0
    natimm = 0
    vacdur = 0
    ref = 0 ## the reference scenario, [0, nsim); if None, plot absolute values
    date_zero = pd.Timestamp('2019-12-31')
    date_startplt = pd.Timestamp('2020-12-15') # pd.Timestamp('2020-12-15') # 
    date_cutoff = pd.Timestamp('2021-06-30') # pd.Timestamp('2021-06-30') # 
    date_startvac = pd.to_datetime(min(startday_1, startday_2), unit='D', origin=date_zero)
    len_hl_a = len(subvac)
    len_subsim = len(subsim)

    if date_startplt > date_zero and (date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        startrow_idx = int((date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) + 1 
    else:
        startrow_idx = 0
    if date_cutoff > date_startplt and (date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        stoprow_idx = int((date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) + 2
    else:
        stoprow_idx = AAA_1.shape[-2]
    
    # %%
    cl_mix = ['darkgreen', 'tab:red', 'xkcd:blue', 'darkturquoise', 'orchid', 'seagreen', 'xkcd:goldenrod',
              'mediumblue', 'darkgreen', 'orangered', 'lime', 'dimgrey']
    ls_mix = ['-', '--']
    style_cycler = cycler(color=np.tile(cl_mix[:len_subsim], len(ls_mix)), linestyle=np.repeat(ls_mix, len_subsim))

    plt.rc('axes', prop_cycle=style_cycler )

    # %%
    ## 6x3 plot
    ## columns: vaccine profiles
    ## rows: low/high supply pairs showing cumulative cases, hospitalization, deaths
    ## colors: vaccine strategies
    ## linestyle: with/without test
    fig, axs = plt.subplots(nrows=6, ncols=len_hl_a, 
                            sharex='col', sharey='row',
                            figsize=(8*len_hl_a, 5.8*6), dpi=120,
                            gridspec_kw={'hspace': 0.10, 'wspace': .15, 
                                         'bottom': 0.08, 'top': 0.95, 'left': 0.13, 'right': 0.68})
    
    if ref is not None:
        A = AAA_1[natimm, norm, vacdur, 0, ref, startrow_idx:stoprow_idx,:]
        base_J = sim.Get_Compartment_All_Ages("J",traj=A)
        base_K = sim.Get_Compartment_All_Ages("K",traj=A)
        base_D = sim.Get_Compartment_All_Ages("D",traj=A) + sim.Get_Compartment_All_Ages("DHOSP",traj=A)
    else:
        A = AAA_1[natimm, norm, vacdur, 0, 0, startrow_idx:stoprow_idx,:]
        base_J = base_K = base_D = 1
    
    ## plot baseline
    for ax in axs.flat:
        if ref is not None:
            ax.plot(A[:,0], [1.0]*A[:,0].shape[0], color='grey', ls='-', label=sim_labels[ref], lw=4.0, alpha=0.5)
    
    ## plot subsim
    for ic, iv in enumerate(subvac):
        v = vac_halflife_a[iv]
        slope = vac_slope_a[iv]
        axs[0,ic].set_title('vaccine halflife: {}\nvaccine slope: {}'.format(v, slope))

        ## with test
        for ia, AAA in enumerate([AAA_1, AAA_3]):
            for iss, s in enumerate(subsim):
                A = AAA[natimm, norm, vacdur, iv, s, startrow_idx:stoprow_idx,:]
                traj_J = sim.Get_Compartment_All_Ages("J",traj=A)
                traj_K = sim.Get_Compartment_All_Ages("K",traj=A)
                traj_D = sim.Get_Compartment_All_Ages("D",traj=A) + sim.Get_Compartment_All_Ages("DHOSP",traj=A)
                print('only Ab-, last shown J, {}: {}'.format(sim_labels[s], traj_J[-1]))
                print('only Ab-, last shown K, {}: {}'.format(sim_labels[s], traj_K[-1]))
                print('only Ab-, last shown D, {}: {}'.format(sim_labels[s], traj_D[-1]))
                axs[0+ia,ic].plot(A[:,0], traj_J / base_J, label='vaccinate only Ab-\n{}'.format(sim_labels[s]), alpha=0.8, lw=4.5)
                axs[2+ia,ic].plot(A[:,0], traj_K / base_K, label='vaccinate only Ab-\n{}'.format(sim_labels[s]), alpha=0.8, lw=4.5)
                axs[4+ia,ic].plot(A[:,0], traj_D / base_D, label='vaccinate only Ab-\n{}'.format(sim_labels[s]), alpha=0.8, lw=4.5)
        ## without test
        for ia, AAA in enumerate([AAA_2, AAA_4]):
            for iss, s in enumerate(subsim):
                A = AAA[natimm, norm, vacdur, iv, s, startrow_idx:stoprow_idx,:]
                traj_J = sim.Get_Compartment_All_Ages("J",traj=A)
                traj_K = sim.Get_Compartment_All_Ages("K",traj=A)
                traj_D = sim.Get_Compartment_All_Ages("D",traj=A) + sim.Get_Compartment_All_Ages("DHOSP",traj=A)
                print('all individuals, last shown J, {}: {}'.format(sim_labels[s], traj_J[-1]))
                print('all individuals, last shown K, {}: {}'.format(sim_labels[s], traj_K[-1]))
                print('all individuals, last shown D, {}: {}'.format(sim_labels[s], traj_D[-1]))
                axs[0+ia,ic].plot(A[:,0], traj_J / base_J, label='vaccinate all individuals\n{}'.format(sim_labels[s]), alpha=0.8, lw=4.5)
                axs[2+ia,ic].plot(A[:,0], traj_K / base_K, label='vaccinate all individuals\n{}'.format(sim_labels[s]), alpha=0.8, lw=4.5)
                axs[4+ia,ic].plot(A[:,0], traj_D / base_D, label='vaccinate all individuals\n{}'.format(sim_labels[s]), alpha=0.8, lw=4.5)
    
    ## legen outside of axes
    axs[0,-1].legend(loc='upper left', bbox_to_anchor=(1.03, 1.0), fontsize=24)

    ## x-axis: date
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
        # ax.axvline(x=(pd.Timestamp('2021-01-01') - date_zero).days, color='black', ls='--', lw=2)
        ax.set_xticks(xdates_ticks)
        ax.set_xticklabels(xdates_labels, fontsize=24)
        ax.set_xlim(left=(date_startplt - date_zero).days, right=(date_cutoff - date_zero).days+1 )
        ax.axvline(x=(date_startvac- date_zero).days, color='grey', ls='--', lw=2)

    ## y-labels
    for ia, ial in enumerate(['LOW SUPPLY', 'AMPLE SUPPLY']):
        if ref is not None:
            axs[0+ia,0].set_ylabel('{}\n\nCumulative cases\nrelative to the\n"{}" strategy'.format(ial, sim_labels[ref]) )
            axs[2+ia,0].set_ylabel('{}\n\nCumulative hospitalizations\nrelative to the\n"{}" strategy'.format(ial, sim_labels[ref]))
            axs[4+ia,0].set_ylabel('{}\n\nCumulative deaths\nrelative to the\n"{}" strategy'.format(ial, sim_labels[ref]))
        else:
            axs[0+ia,0].set_ylabel('{}\n\nCumulative cases'.format(ial))
            axs[2+ia,0].set_ylabel('{}\n\nCumulative hospitalizations'.format(ial) )
            axs[4+ia,0].set_ylabel('{}\n\nCumulative deaths'.format(ial) )

    fig.suptitle('{}, vaccine campaign starts on {}, ends on {} (low supply) or on {} (ample supply), with {} doses per day'.format(loc, 
                        date_startvac.strftime('%d %b %Y'), 
                        pd.to_datetime(endday_1, unit='D', origin=date_zero).strftime('%d %b %Y'),
                        pd.to_datetime(endday_3, unit='D', origin=date_zero).strftime('%d %b %Y'),
                        dpd_1 ), y=0.02)

    if ref is not None: fig.savefig(Path(savefig_dir, '{}-{}_beta-comparetest-i{}_c{}_u{}-ref{}-nsim{}.png'.format(loc, trset, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur], ref, len_subsim )))
    else:               fig.savefig(Path(savefig_dir, '{}-{}_beta-comparetest-i{}_c{}_u{}-abs.png'.format(loc, trset, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur])))






# %%
