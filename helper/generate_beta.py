#!/usr/bin/env python3
# %%
from plot_helper import Set_Axes_Xticks
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from pathlib import Path
from io import BytesIO
from subprocess import run, Popen, PIPE, TimeoutExpired
from matplotlib import rc_file

# plt.rcdefaults()

rc_file(Path(Path.home(), 'Code/covid19-vaccine/helper/matplotlibrc-custom'))


### plot beta over time
def Plot_Betas(beta_list, start_daynum, end_daynum=None, tick_interval=10):
    if end_daynum is None:
        end_daynum = start_daynum + len(beta_list) - 1
    assert (len(beta_list) == end_daynum - start_daynum +1), 'Mismatch length of beta list and daynums.'
    fig, ax = plt.subplots()
    ax = Set_Axes_Xticks(ax, start_daynum, end_daynum, tick_interval=tick_interval)
    ax.plot(range(start_daynum, end_daynum+1), beta_list)
    ax.set_title("Beta")
    ax.grid(b=True, lw=0.8, alpha=0.2)
    plt.tight_layout()
    return fig, ax


### generate betas using sin function
def Generate_Beta_List(partial_beta_list, first_beta_daynum, final_beta_daynum, 
                       mean_beta_start_daynum=None, mean_beta_end_daynum=None, buffer_ndays=30,
                       trough_peak_ratio=0.1, sin_origin_daynum=None,
                       relative_mean_last_betas=1.0):
    # oct-later: A * sin(2 pi * (t – daynum Oct 15) / 365)
    ## | partial_beta_list | 
    ## | changes to reach end-of-buffer value after buffer_ndays | 
    ## | generate using sin function starting from (sin_origin_daynum + 1) |
    len_pbl = len(partial_beta_list)
    assert (len_pbl > 0), 'Beta list must not be empty.'
    if sin_origin_daynum is None: sin_origin_daynum = (pd.Timestamp('2021-10-15') - pd.Timestamp('2019-12-31')).days
    # if sin_origin_daynum is None: sin_origin_daynum = pd.Timestamp('2021-10-15').dayofyear

    if mean_beta_end_daynum is None or mean_beta_start_daynum is None:
        if len_pbl <= 30: mean_beta = np.array(partial_beta_list).mean()
        else: mean_beta = np.array(partial_beta_list[-30:]).mean()
    else:
        assert(mean_beta_start_daynum < mean_beta_end_daynum), 'mean_beta_end_daynum must be higher than mean_beta_start_daynum'
        assert(first_beta_daynum < mean_beta_start_daynum), 'mean_beta_start_daynum must be higher than first_beta_daynum'
        assert(first_beta_daynum < mean_beta_end_daynum), 'mean_beta_end_daynum must be higher than first_beta_daynum'
        if len_pbl <= 30: mean_beta = np.array(partial_beta_list).mean()
        else: mean_beta = np.array(partial_beta_list[(mean_beta_start_daynum - first_beta_daynum):(mean_beta_end_daynum - first_beta_daynum +1)]).mean()

    mean_beta = mean_beta*relative_mean_last_betas
    print(mean_beta)

    final_beta_list = partial_beta_list.copy()
    
    end_buffer_beta = mean_beta if (first_beta_daynum + len_pbl + buffer_ndays - 1) < sin_origin_daynum else (mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*((first_beta_daynum+len_pbl+buffer_ndays) - sin_origin_daynum)/365.0)) )
        
    final_beta_list.extend( list(np.linspace(partial_beta_list[-1], end_buffer_beta, buffer_ndays)) )
    

    if first_beta_daynum + len_pbl + buffer_ndays - 1 < sin_origin_daynum and sin_origin_daynum <  final_beta_daynum:
        final_beta_list.extend([mean_beta]*(sin_origin_daynum - first_beta_daynum - len_pbl - buffer_ndays + 1))

        final_beta_list.extend([mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0)) for t in range(sin_origin_daynum+1, final_beta_daynum+1)])

    else:
        final_beta_list.extend([mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0)) for t in range(first_beta_daynum+len_pbl+buffer_ndays, final_beta_daynum+1)])

    assert(len(final_beta_list) == (final_beta_daynum - first_beta_daynum + 1)), 'len(final_beta_list): {}; it should be {}'.format(len(final_beta_list), final_beta_daynum - first_beta_daynum + 1)
    return final_beta_list

# %%
def Generate_Beta_Fr_Spline(spline_coeff, final_beta_daynum):
    assert (len(spline_coeff) == round((final_beta_daynum-60)/7)), 'Length of spline coefficient list must be round((final_beta_daynum-60)/7).'

    cmd = ["Rscript", "generate_beta_from_spline.R", str(final_beta_daynum)]
    cmd.extend([str(x) for x in spline_coeff])
    ## run command
    proc = run(cmd, capture_output=True, cwd=Path(__file__).parent)
    if proc.returncode == 0:
        beta = np.genfromtxt(BytesIO(proc.stdout), delimiter=" ")
    else: 
        raise Exception("Unable to generate daily betas")
    
    return beta

# %%
if __name__ == "__main__":
    # %%
    out_dir = Path(Path.home(), 'Code/covid19-vaccine/sim_output/beta/RI-20210111/')

    loc = "MA"
    b = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/inference/output/RI/2020_12_15-v6/RI-1215-nonorm_extra_dec06-sampling_100.daily.betas-day-341.csv'))
    p = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/inference/output/RI/2020_12_15-v6/RI-1215-nonorm_extra_dec06-sampling_100.ode.params-day-341.csv'))
    s = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/inference/output/RI/2020_12_15-v6/RI-1215-nonorm_extra_dec06-sampling_100.spline.coeff.betas-day-341.csv'))
    # b = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/inference/output/MA/2020_12_15-v6/MA-1215-nonorm_extra_dec06S-sampling_100.daily.betas-day-341.csv'))
    # p = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/inference/output/MA/2020_12_15-v6/MA-1215-nonorm_extra_dec06S-sampling_100.ode.params-day-341.csv'))
    # s = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/inference/output/MA/2020_12_15-v6/MA-1215-nonorm_extra_dec06S-sampling_100.spline.coeff.betas-day-341.csv'))
    sample_id = None ## set to None to use median
    start_daynum = int(b.columns[0])
    end_daynum = 1096 # start_daynum + 365*2
    
    date_zero = pd.Timestamp('2019-12-31')
    date_start_mean_beta = pd.Timestamp('2020-09-01')
    date_end_mean_beta = pd.Timestamp('2020-11-30')

    # %%
    ## export to space-separated params (i.e. to be used in odesim command)
    if sample_id is not None:
        p_m = p.iloc[sample_id,:].to_dict()
    else:
        p_m = p.median().to_dict()
    # p_m = p.iloc[61,:].to_dict()

    with open(Path(out_dir, '{}-params-{}.txt'.format(loc, 'median' if sample_id is None else 'sampleid_'+str(sample_id))), 'w') as f:
        f.write(str([ '-'+k +' '+ str(round(v,7)) for k, v in p_m.items()]).replace(',', '').replace("'", '')[1:-1])
    # print(*[round(x,7) for x in gen_betas], sep=" ")

    # %%
    ##
    ## base off original betas 
    ##
    rm = 0
    trough_peak_ratio = 0.1
    sin_origin_daynum = (pd.Timestamp('2021-10-15') - date_zero).days
    relative_mean_last_betas = 1.0 # 1.3 # 0.7 # 
    
    if sample_id is not None:
        betas = b.iloc[sample_id,:].to_list() if rm == 0 else b.iloc[sample_id,:(-rm)].to_list()
    else:
        betas = b.median().to_list() if rm == 0 else b.iloc[:,:(-rm)].median().to_list()
    # betas = [x*0.99 for x in betas]

    ## medians only
    gen_betas = Generate_Beta_List(betas, start_daynum, end_daynum, 
                                   trough_peak_ratio=trough_peak_ratio,
                                   mean_beta_start_daynum=(date_start_mean_beta - date_zero).days, 
                                   mean_beta_end_daynum=(date_end_mean_beta - date_zero).days,
                                   relative_mean_last_betas=relative_mean_last_betas,
                                   sin_origin_daynum=sin_origin_daynum,
                                   buffer_ndays=30)
    

    # %%
    ##
    ## plot original betas
    ## 
    fig, ax = plt.subplots()
    ax.plot(range(start_daynum, int(b.columns[-1])+1), b.to_numpy().T, color='grey', lw=1.0, alpha=0.1)

    # ax.plot( np.repeat([int(x) for x in b.columns], b.shape[0]).reshape(-1,b.shape[0]) , b.T, color='grey', lw=1.0, alpha=0.1)
    ax.plot(range(start_daynum, int(b.columns[-1])+1), b.quantile(q=[0.5,0.025,0.975]).to_numpy().T, color='tab:blue')
    ax.axvline(x=int(b.columns[-1])-rm, ls='--', lw=1.5, color='tab:orange', alpha=0.7)
    ax = Set_Axes_Xticks(ax, start_daynum, int(b.columns[-1]) )
    fig.suptitle(loc)
    # ax.set_ylim(top=0.7)
    fig.savefig( Path(out_dir, '{}-beta_start-rm{}.png'.format(loc, rm)) )

    # %%
    ## 
    ## plot generated betas
    ## 
    fig, ax = plt.subplots()

    ## medians only
    ax.plot(range(start_daynum, end_daynum+1), gen_betas, color='tab:blue', lw=2.5, alpha=1.0)

    ## whole posterior
    # ax.plot(range(start_daynum, end_daynum+1), gen_betas.T, color='grey', lw=1.0, alpha=0.1)
    # ax.plot(range(start_daynum, end_daynum+1), np.quantile(gen_betas, q=[0.5,0.025,0.975], axis=0).T, color='tab:blue')

    ax.plot(range(start_daynum, int(b.columns[-1])+1), b.median(), color='tab:green', lw=2.0, alpha=0.6)

    ax.axvline(x=int(b.columns[-1])-rm, ls='--', lw=0.9, color='tab:orange', alpha=0.7)
    ax = Set_Axes_Xticks(ax, start_daynum, end_daynum, tick_interval=90 )
    fig.suptitle('{}, "low season" is {} relative to "mean beta"'.format(loc, relative_mean_last_betas))
    # ax.set_ylim(bottom=0.0, top=0.3)
    fig.tight_layout()
    fig.savefig( Path(out_dir, '{}-med-beta_0p{}-rm{}-{}.png'.format(loc, int(trough_peak_ratio*10), rm, 'median' if sample_id is None else 'sampleid_'+str(sample_id)) ) )


    # %%
    ## export to space-separated list of beta
    with open(Path(out_dir, '{}-med-beta0p{:.0f}-rm{}-{}.txt'.format(loc, trough_peak_ratio*10, rm, 'median' if sample_id is None else 'sampleid_'+str(sample_id))), 'w') as f:
        f.write(str([round(x,7) for x in gen_betas]).replace(',', '')[1:-1])
    # print(*[round(x,7) for x in gen_betas], sep=" ")
    
    
# %%
