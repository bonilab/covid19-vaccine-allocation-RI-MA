#/usr/bin/env python3

import os
import csv
import subprocess
import argparse
# from math import sqrt

## construct command, arguments names and values from dict
def Set_Sim_CMD_From_Dict(cpp_dir,
                output_file= None,
                contact_matrix= False,
                symp_frac_davies= False,
                binary_output= False,
                death_rate_only= False,
                args_dict= None):
    
    ## command as a single string
    # cmd = cpp_dir + "/odesim"  # base command

    # ## construct full command with arguments
    # if output_file is not None: cmd = " ".join([cmd, output_file])
    # else: cmd = " ".join([cmd, "none"])
    # ## contact matrix, available from odesim version 5
    # if contact_matrix: cmd = " ".join([cmd, "-scm"])
    # ## binary output, available from odesim version 5
    # if binary_output: cmd = " ".join([cmd, "-binary-output"])
    # ## death-rate-output, available from odesim version 5
    # if death_rate_only: cmd = " ".join([cmd, "-death-rate-output"]) 

    # for k, v in args_dict.items():
    #     if k == "beta":
    #         cmd = " ".join([cmd, "-beta"] + [str(x) for x in v])
    #     else:
    #         cmd = " ".join([cmd, "-"+k, str(v)])


    ## command as a list for Popen
    cmd = [cpp_dir + "/odesim"]  # base command

    ## construct full command with arguments
    if output_file is not None: cmd.append(output_file)
    else: cmd.append("none")
    ## contact matrix, available from odesim version 5
    if contact_matrix: cmd.append("-scm")
    ## symp-frac-davies, available from odesim version 5
    if symp_frac_davies: cmd.append("-symp-frac-davies")
    ## binary output, available from odesim version 5
    if binary_output: cmd.append("-binary-output")
    ## death-rate-output, available from odesim version 5
    if death_rate_only: cmd.append("-death-rate-output")
    
    ## handle **kwargs
    beta_ls = []
    for k,v in args_dict.items():
        if k == "beta":
            beta_ls.append("-beta")
            beta_ls.extend([str(x) for x in v])
        else:
            # cmd.extend(["-"+k, str(v)])
            cmd.append("-"+k)
            cmd.extend(str(v).split())
    cmd.extend(beta_ls)

    return cmd

##################################################################
## please make sure the simulation is BUILT (i.e. 'make') before running this script
if __name__ == "__main__":
    ## CHANGE this path relative to your HOME directory
    # cpp_dir = os.path.expanduser('/mnt/d/Github/Bio/Projects/covid-19/covid19modeling/cpp-v4-6e-severe-classes/')
    # cpp_dir = os.path.expanduser('~/Code/covid19/cpp-v5-discharges-nonhospdeaths/')
    sim_root_dir = '~/Code/covid19-vaccine/' ## corresponding cpp dir will be appended later

    parser = argparse.ArgumentParser(
        description='Run simulation with betas and parameters from csv files')

    # this is the parser object's way of allowing new arguments to be parsed on the command line
    parser.add_argument(
        '-l', '--loc',
        default='RI',
        type=str,
        help='''Location settings to run the simulation with.
                Values must be either "RI", "MA" or "PA".
                Default: RI''')
    parser.add_argument(
        '-t','--time_final',
        default=0,
        type=int,
        help='''Total days to simulate. If given a non-positive integer, the simulation
                will use the last column header in beta csv as time_final.
                Default: 0''')
    parser.add_argument(
        '-b', '--beta_csv',
        type=str,
        help='''csv file with values for daily betas, each line is for one simulation.
                The first column should be "id", values for betas start from the third column.
                Default: betas.csv''')
    parser.add_argument(
        '-p','--param_csv',
        type=str,
        help='''csv file with values for non-beta values like "dev-hosp-young", "dev-hosp-mid",
                 "dev-hosp-old", "dev-len-hospstay", "dev-icu-frac" each line is for one simulation.
                The first column should be "id", values for these parameters (in that exact order)
                start from the second column.
                Default: params.csv''')
    parser.add_argument(
        '-c','--common_csv',
        type=str,
        help='''csv file with common values that will be used for all simulations.
                The first row should contain headers, i.e. command-line options like "len-nat-immunity",
                "vac1-protect-duration", "vac1-efficacy-halflife", etc.
                The second row should contain values for the corresponding header.
                Default: common.csv''')
    parser.add_argument(
        '--ignore_last_cols',
        type=int,
        default=0,
        help='''number of columns (counted from the end) to ignore from param_csv.
                Default: 0''')
    parser.add_argument(
        '--death_rate_only',
        action='store_true',
        help='''turn on death-rate-output for odesim, available from odesim version 5.
                Default: unset (i.e. False)''')
    parser.add_argument(
        '-s','--sim_version',
        default=6,
        type=int,
        help='''ODE simulation version (possible values: 4,5,6,7).
                Default: 6''')
    parser.add_argument(
        '-m','--steps_per_day',
        default=10,
        type=int,
        help='''integration steps per day, available from odesim version 5.
                Default: 10''')
    parser.add_argument(
        '--davies',
        action='store_true',
        help='''turn on using symp-frac-davies, available from odesim version 5.
                Default: unset (i.e. False)''')
    parser.add_argument(
        '-y','--binary_output',
        action='store_true',
        help='''turn on binary output instead of text, available from odesim version 5.
                Default: unset (i.e. False)''')
    parser.add_argument(
        '-d','--out_dir',
        default='./',
        type=str,
        help='''folder to save output files.
                Default: ./ (i.e. current working directory)''')    
    parser.add_argument(
        'sim_args',
        nargs=argparse.REMAINDER,
        help='''if beta_csv or param_csv is not specified, the remaining arguments will be 
                passed to odesim as-is.''')
    parser.add_argument(
        '--extra_args',
        nargs=argparse.REMAINDER,
        help='''extra arguments to be passed to odesim as-is.''')
    args = parser.parse_args()

    #
    # most common call will be : 
    #         python run_sim.py -b testbetas.csv
    # or : 
    #         python run_sim.py --beta_csv testbetas.csv
    # or :
    #         python run_sim.py -l RI -b testbetas.csv --param_csv testparams.csv
    # or :
    #         python run_sim.py -d output_files/ -l MA -b testbetas.csv -p testparams.csv
    #

    # print(args)
    
    ## sim version
    if args.sim_version == 4: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v4-6e-severe-classes/')
    elif args.sim_version == 5: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v5-discharges-nonhospdeaths/')
    elif args.sim_version == 6: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v6-test-vaccination/')
    elif args.sim_version == 7: cpp_dir = os.path.expanduser(sim_root_dir + 'cpp-v7-prms_by_time/')
    else: raise ValueError('Simulation version must be either 4, 5, or 6.')

    ## check location
    if args.loc not in ["RI", "MA", "PA"] : raise ValueError('Location must be either "RI", "MA" or "PA".')
    
    ## make sure the compiled odesim is up-to-date
    # subprocess.Popen(["make"], cwd=cpp_dir)

    if args.beta_csv is None or args.param_csv is None or args.common_csv is None :
        if len(args.sim_args) < 1: raise Exception('Missing arguments for odesim')
        cmd = ["{}/odesim".format(cpp_dir)]
        cmd.extend( args.sim_args )
        print(cmd)
        # proc = subprocess.Popen(cmd, cwd=os.getcwd(), shell=True) ## cmd as a string
        proc = subprocess.Popen(cmd, cwd=os.getcwd())
        try:
            outs, errs = proc.communicate(timeout=15)
        except subprocess.TimeoutExpired:
            proc.kill()
            outs, errs = proc.communicate()
    else:
        ## the default beta_csv file is betas.csv in the folder where this python script is called
        beta_csv_file = os.path.abspath(args.beta_csv)  # absolute path to beta csv file
        param_csv_file = os.path.abspath(args.param_csv)  # absolute path to param csv file
        common_csv_file = os.path.abspath(args.common_csv)  # absolute path to param csv file
        # print(beta_csv_file)
        # print(param_csv_file)

        ## construct commands and run simulation
        with open(beta_csv_file) as betacsv, open(param_csv_file) as paramcsv, open(common_csv_file) as commoncsv:
            betareader = csv.reader(betacsv, delimiter=',')
            beta_header = next(betareader) # betareader.next()
            paramreader = csv.reader(paramcsv, delimiter=',')
            param_header = next(paramreader) # paramreader.next()
            commonreader = csv.reader(commoncsv, delimiter=',')
            common_header = next(commonreader) # commonreader.next()
            common_values = next(commonreader)

            i = 0
            for line in betareader:

                ## uncomment the three following lines to run ALL combinations of beta and param
                #   paramcsv.seek(0)   ## reset the reader to the begining of the file
                #   next(paramreader) # paramreader.next() ## ignore header
                #   for line_p in paramreader:

                ## comment the following line to run ALL combinations of beta and param
                line_p = next(paramreader) #paramreader.next()
                assert int(line[0]) == int(line_p[0]), 'Mismatch line ID in beta and param csv.'

                # numbetas = len(line)
                tf = args.time_final if args.time_final > 0 else int(beta_header[-1])  # 60+numbetas

                ## construct param dict
                args_dict={'tf': tf, 
                           'loc': args.loc} #,
                        #    'introday': 55}
                
                ## integration steps per day for version 5
                if args.sim_version > 4: args_dict['steps-per-day'] = args.steps_per_day

                ## prepare time-varying params (e.g. tv-hosp-frac-10_1, tv-hosp-frac-10_2, tv-dev-icu-frac-endday)
                tv_dict = dict()
                ## start at column right after id (index), ignore the last 2 columns (reporting rate)
                if args.ignore_last_cols > 0:
                    for k,v in zip(param_header[1:-(args.ignore_last_cols)], line_p[1:-(args.ignore_last_cols)]):
                        if k[:3] == 'tv-': tv_dict[k] = v
                        else: args_dict[k] = v
                else:
                    for k,v in zip(param_header[1:], line_p[1:]):
                        if k[:3] == 'tv-': tv_dict[k] = v
                        else: args_dict[k] = v
                ## process time-varying params
                for k in sorted(tv_dict):
                    k_prefix = k.split(sep='_')[0]
                    if k_prefix in args_dict:
                        args_dict[k_prefix] = args_dict[k_prefix] +' '+ str(tv_dict[k])
                    else:
                        args_dict[k_prefix] = str(tv_dict[k])

                ## add common values to args_dict:
                for k, v in zip(common_header, common_values):
                    args_dict[k] = v
                
                ## add beta list to the argument dict
                args_dict['beta'] = line[1:]

                if args.binary_output: out_file_fmt = '{}/run_output_{}.bin'.format(args.out_dir, line[0])
                else: out_file_fmt = '{}/run_output_{}.tdl'.format(args.out_dir, line[0])

                cmd = Set_Sim_CMD_From_Dict(cpp_dir,
                                            output_file= out_file_fmt,
                                            contact_matrix= args.sim_version > 4,
                                            symp_frac_davies= args.davies,
                                            binary_output= args.binary_output,
                                            death_rate_only= args.death_rate_only,
                                            args_dict= args_dict)
                
                if args.extra_args is not None: 
                    cmd.extend( args.extra_args )
                    # print(cmd)

                if i % 100 == 0:
                    print('run {}'.format(i))

                ## run the simulation
                # subprocess.Popen(cmd, cwd=os.getcwd(), shell=True) ## cmd as a string
                # subprocess.Popen(cmd, cwd=os.getcwd()) ## cmd as a list
                proc = subprocess.Popen(cmd, cwd=os.getcwd())
                try:
                    outs, errs = proc.communicate(timeout=15)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    outs, errs = proc.communicate()

                i += 1
                if i == 1:
                    print(cmd)
                    # break

            print('total runs {}'.format(i))
