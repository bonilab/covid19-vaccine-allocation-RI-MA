#!/usr/bin/env python3

import covid_sim
from numpy import genfromtxt, zeros, fromfile
from pathlib import Path
# from pickle import dump as pcdump
from _pickle import dump as pcdump
import argparse
import os


## read simulation output (tab-delimited text) to sim_all_output ndarray
## the default output_dir is Path(__file__).parent / output-agg
## the default filename format is run_output_{}.tdl, e.g. run_output_0.tdl, run_output_1.tdl .... run_output_<total_sim - 1>.tdl
## the default total_sim is 1000
def Read_txt_To_ndarray(output_dir,
                        total_sim=1000,
                        filename_format="run_output_{}.tdl",
                        traj_ncol=307,
                        binary_output=False,
                        death_rate_only=False,
                        verbose=False):
    # if filename_format is None:
    #     filename_format = "run_{}.txt"
    file_path = os.path.abspath(output_dir + '/' + filename_format.format(0))
    if verbose:
        # print("read 0 :", Path(output_dir, filename_format.format(0)))
        print("read 0 :", file_path )
    
    # if binary_output: tmp_a = fromfile(Path(output_dir, filename_format.format(0))).reshape((-1, traj_ncol))
    # else: tmp_a = genfromtxt(Path(output_dir, filename_format.format(0)), delimiter="\t")
    if binary_output: tmp_a = fromfile( file_path ).reshape((-1, traj_ncol))
    else: tmp_a = genfromtxt(file_path, delimiter="\t")

    if death_rate_only: all_sim_output = zeros((total_sim, tmp_a.shape[0]))
    else: all_sim_output = zeros((total_sim, tmp_a.shape[0], tmp_a.shape[1]))
    all_sim_output[0] = tmp_a

    for i in range(1, total_sim):
        if i % 100 == 0 and verbose:
            print("read", i)

        # if binary_output: all_sim_output[i] = fromfile(Path(output_dir, filename_format.format(i))).reshape((-1, traj_ncol))
        # else: all_sim_output[i] = genfromtxt(Path(output_dir, filename_format.format(i)), delimiter="\t")
        
        file_path = os.path.abspath(output_dir + '/' + filename_format.format(i))
        if binary_output: all_sim_output[i] = fromfile( file_path ).reshape((-1, traj_ncol))
        else: all_sim_output[i] = genfromtxt(file_path, delimiter="\t") 

    if verbose:
        # print("read", i, ":", Path(output_dir, filename_format.format(i)))
        print("read", i, ":", file_path )
    
    print("total files read", i+1)
    
    return all_sim_output

    

###############################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Gather raw simulation output from files and save to a 3D array')

    # this is the parser object's way of allowing new arguments to be parsed on the command line
    parser.add_argument(
        '-d','--dir',
        type=str,
        default=str(Path.cwd()),
        help='''Path the directory with raw output files.
                Default: the current working directory.''')
    parser.add_argument(
        '-t','--total',
        type=int,
        default=1000,
        help='''Total number of raw output files.
                Default: 1000.''')
    parser.add_argument(
        '-f', '--format',
        type=str,
        help='''Filename format. 
                Default: "run_output_{}.tdl" if odesim used text output, "run_output_{}.bin" if binary output
                e.g. run_output_0.tdl, run_output_1.tdl .... run_output_<total_sim - 1>.tdl''')
    parser.add_argument(
        '-s','--sim_version',
        default=6,
        type=int,
        help='''ODE simulation version, either 4, 5, or 6.
                Default: 6''')
    parser.add_argument(
        '-y','--binary_output',
        action='store_true',
        help='''turn on binary output instead of text, available from odesim version 5.
                Default: unset (i.e. False)''')
    parser.add_argument(
        '--death_rate_only',
        action='store_true',
        help='''gather results from odesim runs with death-rate-output option, available from odesim version 5.
                Default: unset (i.e. False)''')
    parser.add_argument(
        '-o', '--output',
        type=str,
        default="all_runs_3darray.dat",
        help='''Name of the output 3D array.
                Default: "all_runs_3darray.dat"''')
    args = parser.parse_args()

    # print(args)
    ## example command
    ## 
    ##      python3.8 gather_output.py -d 20200514/RI/ -t 1000 -f "run_output_{}.tdl" -o 20200514-RI_day151.dat
    ## or
    ##      python3.8 gather_output.py -d 20200514/RI/ -t 1000 --binary_output -o 20200514-RI_day151.dat
    ##

    sim_root_dir = '~/Code/covid19-vaccine/' ## corresponding cpp dir will be appended later

    if args.sim_version == 4: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v4-6e-severe-classes/')
    elif args.sim_version == 5: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v5-discharges-nonhospdeaths/')
    elif args.sim_version == 6: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v6-test-vaccination/')
    else: raise ValueError('Simulation version must be either 4, 5, or 6.')

    if args.format is None: 
        if args.binary_output: file_fmt = "run_output_{}.bin"
        else: file_fmt = "run_output_{}.tdl"
    else:
        file_fmt = args.format

    sim_indices = covid_sim.COVID_SIM(cpp_dir).Get_Indices()
    if args.death_rate_only: traj_ncol = 30
    else: traj_ncol = int(sim_indices['DIMENSION']) + 1 # sim_indices['STARTK'] + sim_indices['NUMAC'] + 1


    ## print( Path(Path.cwd(), args.dir ) )
    ## all_output = Read_txt_To_ndarray(Path(os.getcwd(), args.dir ),
    all_output = Read_txt_To_ndarray(os.path.abspath(os.getcwd() + '/' + args.dir),
                                     total_sim=args.total, 
                                     filename_format=file_fmt,
                                     traj_ncol=traj_ncol,
                                     binary_output=args.binary_output,
                                     death_rate_only=args.death_rate_only)

    ## save aggregated array to file to retrieve later (first run)
    with open(args.output, 'wb') as dump_f:
        pcdump(all_output, dump_f)
    