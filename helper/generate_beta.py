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
    # if first_beta_daynum + len_pbl - 1 < sin_origin_daynum:
    #     final_beta_list.extend([mean_beta]*(sin_origin_daynum - first_beta_daynum - len_pbl + 1))

    #     final_beta_list.extend([mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0)) for t in range(sin_origin_daynum+1, final_beta_daynum+1)])
    #     # final_beta_list.extend([mean_beta + trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0) for t in range(sin_origin_daynum+1, final_beta_daynum+1)])
    # else:
    #     final_beta_list.extend([mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0)) for t in range(first_beta_daynum+len_pbl, final_beta_daynum+1)])
    #     # final_beta_list.extend([mean_beta + trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0) for t in range(first_beta_daynum+len_pbl, final_beta_daynum+1)])

    
    end_buffer_beta = mean_beta if (first_beta_daynum + len_pbl + buffer_ndays - 1) < sin_origin_daynum else (mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*((first_beta_daynum+len_pbl+buffer_ndays) - sin_origin_daynum)/365.0)) )
        
    final_beta_list.extend( list(np.linspace(partial_beta_list[-1], end_buffer_beta, buffer_ndays)) )
    
    # final_beta_list.extend([mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0)) for t in range(first_beta_daynum+len_pbl+buffer_ndays, final_beta_daynum+1)])
    if first_beta_daynum + len_pbl + buffer_ndays - 1 < sin_origin_daynum and sin_origin_daynum <  final_beta_daynum:
        final_beta_list.extend([mean_beta]*(sin_origin_daynum - first_beta_daynum - len_pbl - buffer_ndays + 1))

        final_beta_list.extend([mean_beta + max(0.0, trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0)) for t in range(sin_origin_daynum+1, final_beta_daynum+1)])
        # final_beta_list.extend([mean_beta + trough_peak_ratio*mean_beta*np.sin(2.0*np.pi*(t - sin_origin_daynum)/365.0) for t in range(sin_origin_daynum+1, final_beta_daynum+1)])
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
    out_dir = Path(Path.home(), 'Code/covid19-vaccine/rev_output/beta/20210430-MA_UKBE/')

    loc = "MA"
    b = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/rev_output/beta/20210430-MA_UKBE/MA_UKBE-betas-sample10.csv'), index_col=0)
    p = pd.read_csv(Path(Path.home(), 'Code/covid19-vaccine/rev_output/beta/20210430-MA_UKBE/MA_UKBE-params-sample10.csv'), index_col=0)
    
    sample_id = None ## set to None to use median
    start_daynum = int(b.columns[0])
    end_daynum = start_daynum + 365*2 #365 #1096 # 
    
    date_zero = pd.Timestamp('2019-12-31')
    date_start_mean_beta = pd.Timestamp('2020-09-01')
    date_end_mean_beta = pd.Timestamp('2020-10-31')






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
    relative_mean_last_betas = 0.7 # 1.3 #1.0 # 
    #
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
                                   buffer_ndays=14)
    

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
    fig.savefig( Path(out_dir, '{}-low-beta_0p{}-rm{}-{}.png'.format(loc, int(trough_peak_ratio*10), rm, 'median' if sample_id is None else 'sampleid_'+str(sample_id)) ) )


    # %%
    ## export to space-separated list of beta
    with open(Path(out_dir, '{}-low-beta0p{:.0f}-rm{}-{}.txt'.format(loc, trough_peak_ratio*10, rm, 'median' if sample_id is None else 'sampleid_'+str(sample_id))), 'w') as f:
        f.write(str([round(x,7) for x in gen_betas]).replace(',', '')[1:-1])
    # print(*[round(x,7) for x in gen_betas], sep=" ")
    
    #########################################################
    
    # %%
    tmp = [1.6727976, 1.3091595, 1.0556721, 0.839467, 0.6612355, 0.5318107, 0.4400402, 0.3742544, 0.3336967, 0.2996406, 0.2729291, 0.2547003, 0.2390196, 0.2287049, 0.2211882, 0.2135892, 0.2062044, 0.1985657, 0.1948971, 0.1892483, 0.1783025, 0.1718373, 0.1674553, 0.1618883, 0.1591671, 0.1565791, 0.155222, 0.1547845, 0.1531179, 0.1516463, 0.1469013, 0.1414765, 0.1356835, 0.1281881, 0.1214586, 0.1142545, 0.1079651, 0.0998169, 0.0942869, 0.0902883, 0.0875188, 0.0861613, 0.0851119, 0.0838469, 0.0846169, 0.0855279, 0.0855725, 0.0856781, 0.0851056, 0.0831568, 0.0803501, 0.0763618, 0.0727229, 0.0699567, 0.066197, 0.0623294, 0.0591034, 0.0552461, 0.0524101, 0.048635, 0.0465685, 0.0458828, 0.0468066, 0.0480444, 0.0503944, 0.0528407, 0.0558019, 0.0569346, 0.0581211, 0.0588087, 0.0593474, 0.058193, 0.0566061, 0.056194, 0.0571379, 0.0590063, 0.0616593, 0.0634172, 0.0676003, 0.07469, 0.0857409, 0.0965157, 0.113263, 0.1344348, 0.1585997, 0.1863518, 0.2260575, 0.2611772, 0.2945457, 0.3283191, 0.357635, 0.3868245, 0.4093747, 0.4175377, 0.4246134, 0.4265253, 0.4129889, 0.4171379, 0.4038737, 0.3918221, 0.3705054, 0.3500103, 0.3305635, 0.3021442, 0.293615, 0.2845611, 0.2759113, 0.2699475, 0.2649497, 0.2712738, 0.2789295, 0.2855462, 0.3019285, 0.317321, 0.3335838, 0.3551688, 0.3783127, 0.3904646, 0.4084734, 0.4188742, 0.4410836, 0.4614731, 0.4687238, 0.4800018, 0.4939878, 0.487615, 0.5064311, 0.5060004, 0.5104718, 0.5243046, 0.534116, 0.5293655, 0.5136221, 0.5030991, 0.4887424, 0.4907996, 0.4883806, 0.4959582, 0.4964565, 0.5059589, 0.4983837, 0.49302, 0.4964822, 0.4996463, 0.4999052, 0.4953054, 0.5011131, 0.5030477, 0.5004874, 0.4996323, 0.4951594, 0.4832767, 0.4682172, 0.4563955, 0.4450037, 0.4228513, 0.4111123, 0.395628, 0.3893881, 0.3847026, 0.3812418, 0.3753286, 0.3695844, 0.3793789, 0.383356, 0.3892003, 0.385007, 0.386506, 0.3716624, 0.3765839, 0.3759507, 0.3574669, 0.3559739, 0.3577437, 0.3603862, 0.3612102, 0.3675355, 0.3735665, 0.3779986, 0.3865911, 0.3802716, 0.3724093, 0.3719672, 0.3857572, 0.4022878, 0.4190901, 0.4388562, 0.4524915, 0.4739582, 0.4964955, 0.5164054, 0.5358147, 0.5529294, 0.5704789, 0.5778226, 0.5936818, 0.6159678, 0.6169427, 0.6127916, 0.6235248, 0.617918, 0.6025216, 0.590118, 0.5966963, 0.5890876, 0.584956, 0.5800465, 0.600875, 0.6192757, 0.6417883, 0.6396446, 0.6353462, 0.6359176, 0.6346613, 0.6360578, 0.6331754, 0.6370787, 0.6340793, 0.6227533, 0.6285719, 0.6225454, 0.6215102, 0.6154268, 0.6124175, 0.6120027, 0.6100525, 0.6101937, 0.6151179, 0.6284947, 0.6270344, 0.6247281, 0.623903, 0.6204963, 0.6118679, 0.6098525, 0.6087834, 0.6024977, 0.6051303, 0.6081419, 0.6019591, 0.6034065, 0.5993969, 0.6004755, 0.599124, 0.6028464, 0.6044168, 0.6131865, 0.6183919, 0.6235973, 0.6235772, 0.6256727, 0.6268633, 0.626969, 0.631818, 0.6281357, 0.6239917, 0.6155393, 0.6038295, 0.5981786, 0.5865631, 0.5727228, 0.5607739, 0.5492601, 0.5423473, 0.5444192, 0.54388, 0.5451427, 0.545438, 0.5668876, 0.5975211, 0.6160267, 0.6424269, 0.6639603, 0.6920593, 0.7132619, 0.7379502, 0.7846868, 0.7791358, 0.8053503, 0.8349254, 0.8274932, 0.8274932, 0.8195366, 0.81158, 0.8036233, 0.7956667, 0.7877101, 0.7797535, 0.7717969, 0.7638403, 0.7558837, 0.747927, 0.7399704, 0.7320138, 0.7240572, 0.7161006, 0.708144, 0.7001873, 0.6922307, 0.6842741, 0.6763175, 0.6683609, 0.6604043, 0.6524476, 0.644491, 0.6365344, 0.6285778, 0.6206212, 0.6126646, 0.604708, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5967513, 0.5977786, 0.5988055, 0.5998318, 0.6008571, 0.6018813, 0.6029039, 0.6039248, 0.6049435, 0.6059597, 0.6069733, 0.6079838, 0.608991, 0.6099945, 0.6109942, 0.6119896, 0.6129805, 0.6139666, 0.6149476, 0.6159232, 0.6168931, 0.617857, 0.6188147, 0.6197659, 0.6207102, 0.6216474, 0.6225773, 0.6234995, 0.6244138, 0.6253198, 0.6262175, 0.6271063, 0.6279862, 0.6288569, 0.629718, 0.6305693, 0.6314107, 0.6322417, 0.6330623, 0.6338721, 0.6346709, 0.6354584, 0.6362345, 0.6369989, 0.6377513, 0.6384916, 0.6392195, 0.6399349, 0.6406375, 0.641327, 0.6420033, 0.6426663, 0.6433156, 0.6439511, 0.6445727, 0.64518, 0.6457731, 0.6463516, 0.6469154, 0.6474643, 0.6479982, 0.6485169, 0.6490203, 0.6495082, 0.6499804, 0.6504369, 0.6508775, 0.6513021, 0.6517104, 0.6521025, 0.6524782, 0.6528374, 0.6531799, 0.6535058, 0.6538148, 0.6541069, 0.654382, 0.65464, 0.6548809, 0.6551046, 0.6553109, 0.6554999, 0.6556716, 0.6558257, 0.6559623, 0.6560814, 0.6561829, 0.6562668, 0.6563331, 0.6563817, 0.6564127, 0.6564259, 0.6564215, 0.6563994, 0.6563596, 0.6563022, 0.6562271, 0.6561344, 0.6560241, 0.6558962, 0.6557508, 0.6555879, 0.6554076, 0.6552099, 0.6549949, 0.6547626, 0.6545132, 0.6542466, 0.653963, 0.6536624, 0.6533449, 0.6530107, 0.6526599, 0.6522924, 0.6519085, 0.6515083, 0.6510918, 0.6506592, 0.6502107, 0.6497463, 0.6492662, 0.6487705, 0.6482595, 0.6477331, 0.6471917, 0.6466353, 0.6460641, 0.6454784, 0.6448781, 0.6442637, 0.6436351, 0.6429927, 0.6423365, 0.6416668, 0.6409839, 0.6402878, 0.6395788, 0.6388571]
    plt.plot(tmp)
    ## export to csv
    ## medians only
    #pd.DataFrame(gen_betas).T.to_csv(Path(out_dir, 'med-beta_0p{}.csv'.format(int(trough_peak_ratio*10))), header=np.arange(start_daynum,end_daynum+1), index_label='id' ) 
    
    ## whole posterior
    # pd.DataFrame(gen_betas).to_csv(Path(out_dir, 'beta_0p{}.csv'.format(int(trough_peak_ratio*10))), header=np.arange(start_daynum,end_daynum+1), index_label='id' ) 
    
    







    # %%

    ## whole posterior
    rel_mean_beta = {'med': 1.0, 'low': 0.7, 'high': 1.3}
    
    start_daynum = int(b.columns[0])
    end_daynum = start_daynum + 365*2 #365 #1096 # 
    trs = 'med'
    
    date_zero = pd.Timestamp('2019-12-31')
    date_start_mean_beta = pd.Timestamp('2020-09-01')
    date_end_mean_beta = pd.Timestamp('2020-10-31')

    rm = 0
    trough_peak_ratio = 0.1
    sin_origin_daynum = (pd.Timestamp('2021-10-15') - date_zero).days
    relative_mean_last_betas = rel_mean_beta[trs]
    
    gen_betas = np.array(list(map(lambda betas: Generate_Beta_List(list(betas), start_daynum, end_daynum, 
                                                trough_peak_ratio=trough_peak_ratio, 
                                                mean_beta_start_daynum=(date_start_mean_beta - date_zero).days, 
                                                mean_beta_end_daynum=(date_end_mean_beta - date_zero).days, 
                                                relative_mean_last_betas=relative_mean_last_betas, 
                                                sin_origin_daynum=sin_origin_daynum, 
                                                buffer_ndays=14), 
                                   b.to_numpy() if rm == 0 else b.iloc[:,:(-rm)].to_numpy() ) ))

    # gen_betas = np.array(list(map(lambda betas: Generate_Beta_List(list(betas), start_daynum, end_daynum, trough_peak_ratio=trough_peak_ratio), b.iloc[:,:(-rm)].to_numpy() ) ))

    tmp = pd.DataFrame(gen_betas, index=b.index, columns=np.arange(start_daynum, end_daynum+1))

    tmp.to_csv(Path(out_dir, 'MA_UKBE-gen_betas-sample10-{}.csv'.format(trs)))
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
    ax.plot(range(start_daynum, end_daynum+1), gen_betas.T, color='grey', lw=1.0, alpha=0.1)

    ax.plot(range(start_daynum, end_daynum+1), np.quantile(gen_betas, q=[0.5,0.025,0.975], axis=0).T, color='tab:blue')
    ax.axvline(x=int(b.columns[-1])-rm, ls='--', lw=1.5, color='tab:orange', alpha=0.7)
    ax = Set_Axes_Xticks(ax, start_daynum, int(b.columns[-1]) )
    fig.suptitle(loc)
    # ax.set_ylim(top=0.7)
    fig.savefig( Path(out_dir, '{}-beta_{}-rm{}.png'.format(loc, trs, rm)) )

    # %%
    ##
    ## using rolling mean
    ##
    rm = 0
    trough_peak_ratio = 0.1
    rel_b = 0.82 # 0.95 # 0.70 #  ## overall new betas should be rel_b of original betas
    len_buffer = 60 ## days of using rolling mean to generate new betas right after the original 
    b_m_ls = b.iloc[:,:(-rm)].median().to_list()
    # b_m_ls.extend( b.median().multiply(rel_b).rolling(3).mean().to_list()[-1:-(len_buffer+1):-1] )
    b_m_ls.extend( b.iloc[:,:(-rm)].median().multiply(rel_b).rolling(3).mean().to_list()[-1:-(len_buffer+1):-1] )
    # b_m_ls = b.iloc[70,:(-rm)].to_list()
    # b_m_ls.extend( b.iloc[70,:].multiply(rel_b).rolling(3).mean().to_list()[-1:-(len_buffer+1):-1] )

    ## medians only
    gen_betas = Generate_Beta_List(b_m_ls, start_daynum, end_daynum, trough_peak_ratio=trough_peak_ratio,
                                    relative_mean_last_betas=1.0)

      

    



    # %%    
    fig, ax = Plot_Betas(betas, 61, tick_interval=90)
    plt.tight_layout()
    # fig.savefig( Path(out_dir, 'beta_start.png') )

    # %%
    ## using spline
    ##
    len_spl = round((end_daynum-60)/7)
    sp_m_ls = s.median().to_list()
    # sp_m_ls.extend( s.median().rolling(4).mean().to_list()[-1:(-len(sp_m_ls)+16):-1]  )
    sp_m_ls.extend( s.median().rolling(4).mean().to_list()[-1:-17:-1]  )

    beta_fr_spl = Generate_Beta_Fr_Spline(sp_m_ls, len(sp_m_ls)*7+60)

    fig, ax = Plot_Betas(list(beta_fr_spl), 61, tick_interval=90)
    ax.plot(range(61,342), b.median())
    # ax.set_ylim(bottom=0.0, top=1.0)
    plt.tight_layout()

    # fig.savefig( Path(out_dir, 'beta0p{:.0f}.png'.format(trough_peak_ratio*10)) )

    # %%
    ##
    ## using rolling mean
    ##
    rm = 7
    trough_peak_ratio = 0.1
    rel_b = 0.90  ## overall new betas should be rel_b of original betas
    len_buffer = 60 ## days of using rolling mean to generate new betas right after the original 
    b_m_ls = b.iloc[:,:(-rm)].median().to_list()
    b_m_ls.extend( b.median().multiply(rel_b).rolling(7).mean().to_list()[-1:-(len_buffer+1):-1] )

    ## medians only
    b_ls = Generate_Beta_List(b_m_ls, start_daynum, end_daynum, trough_peak_ratio=trough_peak_ratio,
                                    relative_mean_last_betas=1.0)

    # b_ls_0 = Generate_Beta_List(b.iloc[:,:(-rm)].median().to_list(), start_daynum, end_daynum, 
    #                           trough_peak_ratio=trough_peak_ratio, relative_mean_last_betas=1.0,
    #                           mean_beta_start_daynum=pd.Timestamp('2020-07-01').dayofyear,
    #                           mean_beta_end_daynum=pd.Timestamp('2020-07-31').dayofyear )
    
    b_ls_0 = Generate_Beta_List(b_m_ls, start_daynum, end_daynum, 
                              trough_peak_ratio=trough_peak_ratio, relative_mean_last_betas=1.0,
                              mean_beta_start_daynum=int(b.columns[-(len_buffer+rm+31)]),
                              mean_beta_end_daynum=int(b.columns[-(len_buffer+rm+1)]) )

    fig, ax = Plot_Betas(b_ls, 61, tick_interval=90)
    # ax.plot(range(61,61+len(b_ls_0)), b_ls_0, alpha=0.9 )
    # ax.plot(range(61,61+len(b_ls_0)), pd.Series(b_ls_0).rolling(7).mean() )
    ax.plot(range(61,61+len(b_m_ls)), b_m_ls, alpha=0.8 )
    ax.plot(range(61,342), b.median(), alpha=0.5)
    # ax.set_ylim(bottom=0.0, top=1.0)
    plt.tight_layout()
    # %%
    rm = 7
    trough_peak_ratio = 0.1
    rel_b = 0.90  ## overall new betas should be rel_b of original betas
    len_buffer = 60 ## days of using rolling mean to generate new betas right after the original 
    b_m_ls = b.iloc[:,:(-rm)].median().to_list()
    b_m_ls.extend( b.median().multiply(rel_b).rolling(7).mean().to_list()[-1:-(len_buffer+1):-1] )

    ## medians only
    b_ls = Generate_Beta_List(b_m_ls, start_daynum, end_daynum, trough_peak_ratio=trough_peak_ratio,
                                    relative_mean_last_betas=1.0)

    
    fig, ax = Plot_Betas(b_ls, 61, tick_interval=90)
    
    plt.tight_layout()
# %%
