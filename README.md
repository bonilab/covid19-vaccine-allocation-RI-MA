# covid19-vaccine-allocation-RI-MA
Model evaluation of SARS-CoV-2 vaccine allocation in Rhode Island and Massachusets in Spring 2021

Preprint: [Optimal SARS-CoV-2 vaccine allocation using real-time seroprevalence estimates in Rhode Island and Massachusetts](https://doi.org/10.1101/2021.01.12.21249694) 

## 1. Compile C++ code
```shell
cd cpp-v6-test-vaccination
make
```

## 2. Reproduce trajectories
Commands necessary for reproducing trajectories under two vaccine supplies and three transmission settings scenarios for Rhode Island (RI) and Massachusetts (MA) are available in `sim_output`. The easiest way to reproduce trajectories shown in Figure 3, 5-7 (and similar Supplementary Figures) is to modify and run `RI-tot_gen_cmd_and_run.sh` for RI (and `MA-tot_gen_cmd_and_run.sh` for MA).

```shell
cd sim_output
chmod +x RI-tot_gen_cmd_and_run.sh
./RI-tot_gen_cmd_and_run.sh
```

Scripts to reproduce trajectories used to make violin plots in Figure 4 and Figure S13 are available in `20201226-med-inexgrp-RI-tot_300k` and `20210109-med-inexgrp-MA-tot_1800k`, respectively. Usage:

```shell
## for RI
cd 20201226-med-inexgrp-RI-tot_300k
chmod +x run-in_ex_group-ratio-full_doses.sh
mkdir s370_e429_d5000 
./run-in_ex_group-ratio-full_doses.sh -o s370_e429_d5000
```

```shell
## for MA
cd 20210109-med-inexgrp-MA-tot_1800k
chmod +x run-in_ex_group-ratio-full_doses.sh
mkdir s370_e429_d30000 
./run-in_ex_group-ratio-full_doses.sh -o s370_e429_d30000
```

**Note:** raw trajectories files are available on request.

## 3. Reproduce figures
Python3 scripts needed to reproduce main text figures (and similar for those in the supplement) are available in `helper` folder:

- Figure 2: `plot_efficacy.py`
- Figure 3: `plot_low_med_high_transmission.py`
- Figure 4: `plot_vac_in_ex_group.py`
- Figure 5: `plot_vac_all_combine.py`
- Figure 6: `plot_vac_endpoint_withZ_combine.py`
- Figure 7: `plot_vac_serotest.py`
