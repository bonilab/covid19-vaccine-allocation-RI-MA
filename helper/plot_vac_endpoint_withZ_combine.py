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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
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

    loc = 'MA'
    trset = 'high'
    savefig_dir = Path(cpp_dir.parent, 'sim_output/20210105-figures')
    out_dir_1 = Path(cpp_dir.parent, 'sim_output/20201228-{}-sens204970-{}-tot_300k'.format(trset, loc))
    out_dir_2 = Path(cpp_dir.parent, 'sim_output/20210105-{}-sens204970-{}-tot_1800k'.format(trset, loc))
    ## vaccination stratgies:
    # sim_labels = ['no vaccine', 'random', 
    #               'random 10/90 among 20-49 and 70+',
    #               'random 20/80 among 20-49 and 70+',
    #               'random 30/70 among 20-49 and 70+',
    #               'random 40/60 among 20-49 and 70+',
    #               'random 50/50 among 20-49 and 70+',
    #               'random 60/40 among 20-49 and 70+',
    #               'random 70/30 among 20-49 and 70+',
    #               'random 80/20 among 20-49 and 70+',
    #               'random 90/10 among 20-49 and 70+']
    sim_labels = ['no vaccine', 'random',  
                  '10/90',
                  '20/80',
                  '30/70',
                  '40/60',
                  '50/50',
                  '60/40',
                  '70/30',
                  '80/20',
                  '90/10']
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
    # sim_labels = ['no vaccine', 'random', 
    #               'random 10/90 among 20-49 and 70+',
    #               'random 20/80 among 20-49 and 70+',
    #               'random 30/70 among 20-49 and 70+',
    #               'random 40/60 among 20-49 and 70+',
    #               'random 50/50 among 20-49 and 70+',
    #               'random 60/40 among 20-49 and 70+',
    #               'random 70/30 among 20-49 and 70+',
    #               'random 80/20 among 20-49 and 70+',
    #               'random 90/10 among 20-49 and 70+']
    # nat_imm_a=[540]
    # normalcyday_a=[1098]
    # vac_duration_a=[540]
    # vac_halflife_a=[180, 180, 180, 360, 360, 540, 540]
    # vac_slope_a =  [2,   3,   4,   1.5, 2,   1.5, 2]
    # %%
    subsim = [x for x in range(2,nsim)][::-1] # [1] + [x for x in range(3,nsim)] #
    subvac = [2, 5] # [x for x in range(0, len(vac_halflife_a)) ]
    len_subsim = len(subsim)
    norm = 0
    natimm = 0
    vacdur = 0
    ref = 1 ## the reference scenario, [0, nsim); if None, plot absolute values
    date_zero = pd.Timestamp('2019-12-31')
    date_startplt = pd.Timestamp('2020-12-15') # pd.Timestamp('2020-03-01') # 
    date_cutoff = pd.Timestamp('2021-06-30') # pd.Timestamp('2022-12-31') # 
    date_startvac = pd.to_datetime(min(startday_1, startday_2), unit='D', origin=date_zero)
    len_hl_a = len(subvac)

    if date_startplt > date_zero and (date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        startrow_idx = int((date_startplt - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) 
    else:
        startrow_idx = 0
    if date_cutoff > date_startplt and (date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] < AAA_1.shape[-2]:
        stoprow_idx = int((date_cutoff - date_zero).days - AAA_1[0,0,0,0, 0,0,0] ) + 1
    else:
        stoprow_idx = AAA_1.shape[-2]
    
    
    # %%
    ## include cumulative vaccinations
    ##
    ### colors
    len_subsim = len(subsim)

    newcmap = cm.get_cmap('Spectral', len_subsim*2).reversed()
    cl_subsim = [newcmap(0)] + [newcmap(x/len_subsim) for x in range(1,len_subsim-1)] + [newcmap(1)]

    # cl_subsim[len_subsim//2] = newcmap(x/len_subsim)
    cl_mix = np.repeat(cl_subsim, 2, axis=0)
    ls_mix = ['--', '-'] * len_subsim
    

    style_cycler = cycler(color=cl_mix, linestyle=ls_mix)
    # style_cycler = cycler(color=cl_t10_css[:nsim])

    plt.rc('axes', prop_cycle=style_cycler )

    
    ## 
    ## cumulativehospitalization, deaths  
    ##
    ndoses_1 = (endday_1 - startday_1 + 1)*dpd_1
    ndoses_2 = (endday_2 - startday_2 + 1)*dpd_2
    fig, axs = plt.subplots(nrows=2, ncols=len_hl_a*2, 
                            sharex=False, sharey='col',
                            figsize=(11*len_hl_a*2, 25), dpi=120,
                            gridspec_kw={'hspace': 0.24, 'wspace': 0.35, 
                                         'bottom': 0.12, 'top': 0.94, 'left': 0.05, 'right': 0.81,
                                         'width_ratios': [1,2.4,1,2.4]})

    ymax = ymin = 0
    for ir, AAA, ndoses in zip(range(2), [AAA_1, AAA_2], [ndoses_1, ndoses_2]):
        for ic, iv in enumerate(subvac):
            v = vac_halflife_a[iv]
            xlab = []
            
            ## first row: percentages of hospitalization, deaths
            ## AAA_1
            if ref is not None:
                base_J = sim.Get_Compartment_All_Ages("J",traj=AAA[natimm, norm, vacdur, iv, ref, startrow_idx:stoprow_idx,:])
                base_K = sim.Get_Compartment_All_Ages("K",traj=AAA[natimm, norm, vacdur, iv, ref, startrow_idx:stoprow_idx,:])
                base_D = sim.Get_Compartment_All_Ages("D",traj=AAA[natimm, norm, vacdur, iv, ref, startrow_idx:stoprow_idx,:]) + \
                        sim.Get_Compartment_All_Ages("DHOSP",traj=AAA[natimm, norm, vacdur,iv, ref, startrow_idx:stoprow_idx,:])
            else:
                base_J = base_K = base_D = 1
            for iss, s in enumerate(subsim):
                A = AAA[natimm, norm, vacdur, iv, s,startrow_idx:stoprow_idx,:]
                tot_Z_20_49 = sim.Get_Compartment_Age_Group("JZ_1", 2, traj=A) + sim.Get_Compartment_Age_Group("JZ_1", 3, traj=A) + sim.Get_Compartment_Age_Group("JZ_1", 4, traj=A)
                tot_Z_70_80 = sim.Get_Compartment_Age_Group("JZ_1", 7, traj=A) + sim.Get_Compartment_Age_Group("JZ_1", 8, traj=A) 
                reduction_J = 100.0*( ((sim.Get_Compartment_All_Ages("J",traj=A)) / base_J)[-1] - 1.0)
                reduction_K = 100.0*( ((sim.Get_Compartment_All_Ages("K",traj=A)) / base_K)[-1] - 1.0)
                reduction_D = 100.0*( ((sim.Get_Compartment_All_Ages("D",traj=A) + sim.Get_Compartment_All_Ages("DHOSP",traj=A)) / base_D)[-1] - 1.0)
                # ymax = max(ymax, max(reduction_D, reduction_K, reduction_J))
                # ymin = min(ymin, min(reduction_D, reduction_K, reduction_J))
                ymax = max(ymax, max(reduction_D, reduction_K))
                ymin = min(ymin, min(reduction_D, reduction_K))
                axs[ir,ic*2].plot(A[:,0], tot_Z_20_49, lw=5, alpha=0.90, label='{} - 20-49 group'.format(sim_labels[s]))
                axs[ir,ic*2].plot(A[:,0], tot_Z_70_80, lw=5, alpha=0.90, label='{} - 70+ group'.format(sim_labels[s]))
                if iss == 0:
                    # axs[ir,ic*2+1].plot(len_subsim - 1 - iss, reduction_J, ls='', marker='s', fillstyle='none', ms=20, mew=3, color='black', label='cumulative symptomatic cases')
                    axs[ir,ic*2+1].plot(len_subsim - 1 - iss, reduction_K, ls='', marker='o', fillstyle='none', ms=20, mew=3, color='black', label='cumulative hospitalizations')
                    axs[ir,ic*2+1].plot(len_subsim - 1 - iss, reduction_D, ls='', marker='x', fillstyle='none', ms=20, mew=3, color='black', label='cumulative deaths')
                else:
                    # axs[ir,ic*2+1].plot(len_subsim - 1 - iss, reduction_J, ls='', marker='s', fillstyle='none', ms=20, mew=3, color='black')
                    axs[ir,ic*2+1].plot(len_subsim - 1 - iss, reduction_K, ls='', marker='o', fillstyle='none', ms=20, mew=3, color='black')
                    axs[ir,ic*2+1].plot(len_subsim - 1 - iss, reduction_D, ls='', marker='x', fillstyle='none', ms=20, mew=3, color='black')
                xlab.append('{:.1%}\n{:.1%}'.format(round(tot_Z_20_49[-1])/ndoses, round(tot_Z_70_80[-1])/ndoses))
            axs[ir,ic*2+1].set_xticks(range(len_subsim))
            axs[ir,ic*2+1].set_xticklabels(xlab[::-1])
            axs[ir,ic*2+1].axhline(y=0.0, color='grey', ls='--', lw=1.2)
            if ir == 0:
                axs[ir,ic*2].set_title('vaccine halflife: {}\nvaccine slope: {}'.format(v, vac_slope_a[iv]))
                axs[ir,ic*2+1].set_title('vaccine halflife: {}\nvaccine slope: {}'.format(v, vac_slope_a[iv]))

    # ymin = min(int(ymin), -5)
    # ymax = max(int(ymax), 3)
    ymin = min(int(ymin), -5)
    ymax = max(int(ymax), 2)
    if ref is not None:
        for ax in axs[:,[1,-1]].flat:
            ax.set_yticks(range(ymin, ymax, 1))
            ax.set_yticklabels(['{}%'.format(y) for y in range(ymin, ymax, 1)])
            ax.set_ylabel('Percentage compared to the "{}" strategy\non {}'.format(sim_labels[ref], date_cutoff.strftime('%d %b %Y')) )
    else:
        for ax in axs[:,[1,-1]].flat:
            ax.set_ylabel('As of {}'.format(date_cutoff.strftime('%d %b %Y')) )
    
    for ax in axs[:,[1,-1]].flat:
        for iss, txt in zip(range(len_subsim), ax.get_xticklabels()):
            txt.set_color(cl_subsim[len_subsim - 1 - iss])
            txt.set_fontweight('bold')
        # ax.text(-0.6, ymin-2.71, '20-49 group: \n70+ group: ', {'ha': 'right'}, fontsize=24)
        ax.text(-0.6, ymin-2.67, '20-49 group: \n70+ group: ', {'ha': 'right'}, fontsize=24)
        ax.set_xlabel('Percentage of total doses received by each group\n by the end of the campaign', labelpad=22 )
    
    for ax in axs[0,[1,-1]].flat:
        ax.text(len(subsim) - 0.7, ymin-0.1, '10-day campaign\ntotal {} doses'.format(ndoses_1), {'ha': 'right'}, fontsize=22)
    for ax in axs[1,[1,-1]].flat:
        ax.text(len(subsim) - 0.7, ymin-0.1, '60-day campaign\ntotal {} doses'.format(ndoses_2), {'ha': 'right'}, fontsize=22)
    
    for ax in axs[:,[0,2]].flat:
        ax = Set_Axes_Xticks(ax, (date_startvac - date_zero).days, (date_cutoff- date_zero).days, tick_interval=20)
        ax.set_ylabel('Cumulative vaccinations', fontsize=24)
        ax.set_xlim(right=(date_cutoff- date_zero).days - 100)
    
    axs[0,-1].legend(loc='upper left', bbox_to_anchor=(1.03, 1.0), fontsize=25)
    axs[0,-2].legend(loc='upper left', bbox_to_anchor=(4.08, 0.84), fontsize=25)

    fig.suptitle('{}, vaccine campaign starts on {}, ends on {} (low supply) or on {} (ample supply), with {} doses per day'.format(loc, 
                        date_startvac.strftime('%d %b %Y'), 
                        pd.to_datetime(endday_1, unit='D', origin=date_zero).strftime('%d %b %Y'),
                        pd.to_datetime(endday_2, unit='D', origin=date_zero).strftime('%d %b %Y'),
                        dpd_1 ), y=0.02)
            
    # if ref is not None: fig.savefig(Path(savefig_dir, '{}-med_beta-rand20_49_70-i{}_c{}_u{}-combine_withZJ-ref{}.png'.format(loc, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur], ref )))
    # else:               fig.savefig(Path(savefig_dir, '{}-med_beta-rand20_49_70-i{}_c{}_u{}-combine-_withZJ-abs.png'.format(loc, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur])))

    if ref is not None: fig.savefig(Path(savefig_dir, '{}-{}_beta-rand20_49_70-i{}_c{}_u{}-combine_withZ-ref{}.png'.format(loc, trset, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur], ref )))
    else:               fig.savefig(Path(savefig_dir, '{}-{}_beta-rand20_49_70-i{}_c{}_u{}-combine-_withZ-abs.png'.format(loc, trset, nat_imm_a[natimm], normalcyday_a[norm], vac_duration_a[vacdur])))










# %%

