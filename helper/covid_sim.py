#!/usr/bin/env python3

# %%
from numpy import fromfile, frombuffer, genfromtxt, sum as npsum
from io import BytesIO
from subprocess import run, Popen, PIPE, TimeoutExpired
# from shlex import split as shsplit
from pathlib import Path  ## python3.8
# from os import path, getcwd ## python2.7

# %%
class COVID_SIM:
    def __init__(self, cpp_directory):
        self.cpp_dir = Path(cpp_directory).expanduser()
        ## tuple showing the validity lifetime of params
        ## e.g. (4,5) means the param is available only in version 4
        ## last updated: 2020-06-01 (as of odesim version 5)
        self.param_sim_version =   {"introday": (4,7),
                                    "tf": (4,7),
                                    "symp-frac": (4,5),
                                    "dev-hosp-young": (4,5),
                                    "dev-hosp-mid": (4,5),
                                    "dev-hosp-old": (4,5),
                                    "dev-icu-frac": (4,7),
                                    "dev-icu-frac-phase2": (5,7),
                                    "dev-icu-frac-phase2beginday": (5,7),
                                    "dev-len-hospstay": (4,7),
                                    "dev-ventdeath-mid": (4,7),
                                    "rel-beta-hosp": (4,7),
                                    "beta": (4,7),
                                    "death-prob-home-60": (5,7),
                                    "death-prob-home-70": (5,7),
                                    "death-prob-home-80": (5,7),
                                    "death-prob-nonicu-80": (5,7),
                                    "earlymarch-hosp-factor": (5,7),
                                    "hosp-frac-10": (5,7),
                                    "hosp-frac-20": (5,7),
                                    "hosp-frac-30": (5,7),
                                    "hosp-frac-40": (5,7),
                                    "hosp-frac-50": (5,7),
                                    "hosp-frac-60": (5,7),
                                    "hosp-frac-70": (5,7),
                                    "hosp-frac-80": (5,7),
                                    "mean-time-vent": (5,7),
                                    "prob-icu-vent": (5,7),
                                    "self-isolation-factor": (5,7),
                                    "susc-0-20": (5,7),
                                    "susc-60-100": (5,7),
                                    "time-symp-to-hosp": (5,7),
                                    "min-mixinglevel-00": (5,7),
                                    "min-mixinglevel-10": (5,7),
                                    "min-mixinglevel-20": (5,7),
                                    "min-mixinglevel-30": (5,7),
                                    "min-mixinglevel-40": (5,7),
                                    "min-mixinglevel-50": (5,7),
                                    "min-mixinglevel-60": (5,7),
                                    "min-mixinglevel-70": (5,7),
                                    "min-mixinglevel-80": (5,7),
                                    "contact-rate-10": (5,6),
                                    "contact-rate-20": (5,6),
                                    "contact-rate-30": (5,6),
                                    "contact-rate-40": (5,6),
                                    "contact-rate-50": (5,6),
                                    "contact-rate-60": (5,6),
                                    "contact-rate-70": (5,6),
                                    "contact-rate-80": (5,6),
                                    "firstlockdown-endday": (5,7),
                                    "contact-rate-postld-10": (5,6),
                                    "contact-rate-postld-20": (5,6),
                                    "contact-rate-postld-30": (5,6),
                                    "contact-rate-postld-40": (5,6),
                                    "contact-rate-postld-50": (5,6),
                                    "contact-rate-postld-60": (5,6),
                                    "contact-rate-postld-70": (5,6),
                                    "contact-rate-postld-80": (5,6),
                                    "contact-coeff-00": (6,7),
                                    "contact-coeff-10": (6,7),
                                    "contact-coeff-20": (6,7),
                                    "contact-coeff-30": (6,7),
                                    "contact-coeff-40": (6,7),
                                    "contact-coeff-50": (6,7),
                                    "contact-coeff-60": (6,7),
                                    "contact-coeff-70": (6,7),
                                    "contact-coeff-80": (6,7),
                                    "contact-coeff-postld-00": (6,7),
                                    "contact-coeff-postld-10": (6,7),
                                    "contact-coeff-postld-20": (6,7),
                                    "contact-coeff-postld-30": (6,7),
                                    "contact-coeff-postld-40": (6,7),
                                    "contact-coeff-postld-50": (6,7),
                                    "contact-coeff-postld-60": (6,7),
                                    "contact-coeff-postld-70": (6,7),
                                    "contact-coeff-postld-80": (6,7),
                                    "normalcy-beginday": (6,7),
                                    "contact-coeff-norm-00": (6,7),
                                    "contact-coeff-norm-10": (6,7),
                                    "contact-coeff-norm-20": (6,7),
                                    "contact-coeff-norm-30": (6,7),
                                    "contact-coeff-norm-40": (6,7),
                                    "contact-coeff-norm-50": (6,7),
                                    "contact-coeff-norm-60": (6,7),
                                    "contact-coeff-norm-70": (6,7),
                                    "contact-coeff-norm-80": (6,7),
                                    "len-nat-immunity": (6,7),
                                    "vac1-protect-duration": (6,7),
                                    "vac1-efficacy-halflife": (6,7),
                                    "vac1-efficacy-slope": (6,7),
                                    "vac1-phase1-beginday": (6,7),
                                    "vac1-phase1-endday": (6,7),
                                    "vac1-phase1-dpd": (6,7),
                                    "vac1-phase2-beginday": (6,7),
                                    "vac1-phase2-endday": (6,7),
                                    "vac1-frac-10": (6,7),
                                    "vac1-frac-20": (6,7),
                                    "vac1-frac-30": (6,7),
                                    "vac1-frac-40": (6,7),
                                    "vac1-frac-50": (6,7),
                                    "vac1-frac-60": (6,7),
                                    "vac1-frac-70": (6,7),
                                    "vac1-frac-80": (6,7),
                                    "vac2-protect-duration": (6,7),
                                    "vac2-efficacy-halflife": (6,7),
                                    "vac2-efficacy-slope": (6,7),
                                    "vac2-phase1-beginday": (6,7),
                                    "vac2-phase1-endday": (6,7),
                                    "vac2-phase1-dpd": (6,7),
                                    "vac2-phase2-beginday": (6,7),
                                    "vac2-phase2-endday": (6,7),
                                    "vac2-frac-10": (6,7),
                                    "vac2-frac-20": (6,7),
                                    "vac2-frac-30": (6,7),
                                    "vac2-frac-40": (6,7),
                                    "vac2-frac-50": (6,7),
                                    "vac2-frac-60": (6,7),
                                    "vac2-frac-70": (6,7),
                                    "vac2-frac-80": (6,7)
                                }
        ## base unset param_sim 
        self.param_sim =   {"introday": None,
                            "tf": None,
                            "symp-frac": None,
                            "dev-hosp-young": None,
                            "dev-hosp-mid": None,
                            "dev-hosp-old": None,
                            "dev-icu-frac": None,
                            "dev-icu-frac-phase2": None,
                            "dev-icu-frac-phase2beginday": None,
                            "dev-len-hospstay": None,
                            "dev-ventdeath-mid": None,
                            "rel-beta-hosp": None,
                            "beta": None,
                            "death-prob-home-60": None,
                            "death-prob-home-70": None,
                            "death-prob-home-80": None,
                            "death-prob-home-80": None,
                            "earlymarch-hosp-factor": None,
                            "hosp-frac-10": None,
                            "hosp-frac-20": None,
                            "hosp-frac-30": None,
                            "hosp-frac-40": None,
                            "hosp-frac-50": None,
                            "hosp-frac-60": None,
                            "hosp-frac-70": None,
                            "hosp-frac-80": None,
                            "mean-time-vent": None,
                            "prob-icu-vent": None,
                            "self-isolation-factor": None,
                            "susc-0-20": None,
                            "susc-60-100": None,
                            "time-symp-to-hosp": None,
                            "min-mixinglevel-00": None,
                            "min-mixinglevel-10": None,
                            "min-mixinglevel-20": None,
                            "min-mixinglevel-30": None,
                            "min-mixinglevel-40": None,
                            "min-mixinglevel-50": None,
                            "min-mixinglevel-60": None,
                            "min-mixinglevel-70": None,
                            "min-mixinglevel-80": None,
                            "contact-rate-10": None,
                            "contact-rate-20": None,
                            "contact-rate-30": None,
                            "contact-rate-40": None,
                            "contact-rate-50": None,
                            "contact-rate-60": None,
                            "contact-rate-70": None,
                            "contact-rate-80": None,
                            "firstlockdown-endday": None,
                            "contact-rate-postld-10": None,
                            "contact-rate-postld-20": None,
                            "contact-rate-postld-30": None,
                            "contact-rate-postld-40": None,
                            "contact-rate-postld-50": None,
                            "contact-rate-postld-60": None,
                            "contact-rate-postld-70": None,
                            "contact-rate-postld-80": None,
                            "contact-coeff-00": None,
                            "contact-coeff-10": None,
                            "contact-coeff-20": None,
                            "contact-coeff-30": None,
                            "contact-coeff-40": None,
                            "contact-coeff-50": None,
                            "contact-coeff-60": None,
                            "contact-coeff-70": None,
                            "contact-coeff-80": None,
                            "contact-coeff-postld-00": None,
                            "contact-coeff-postld-10": None,
                            "contact-coeff-postld-20": None,
                            "contact-coeff-postld-30": None,
                            "contact-coeff-postld-40": None,
                            "contact-coeff-postld-50": None,
                            "contact-coeff-postld-60": None,
                            "contact-coeff-postld-70": None,
                            "contact-coeff-postld-80": None,
                            "normalcy-beginday": None,
                            "contact-coeff-norm-00": None,
                            "contact-coeff-norm-10": None,
                            "contact-coeff-norm-20": None,
                            "contact-coeff-norm-30": None,
                            "contact-coeff-norm-40": None,
                            "contact-coeff-norm-50": None,
                            "contact-coeff-norm-60": None,
                            "contact-coeff-norm-70": None,
                            "contact-coeff-norm-80": None,
                            "len-nat-immunity": None,
                            "vac1-protect-duration": None,
                            "vac1-efficacy-halflife": None,
                            "vac1-efficacy-slope": None,
                            "vac1-phase1-beginday": None,
                            "vac1-phase1-endday": None,
                            "vac1-phase1-dpd": None,
                            "vac1-phase2-beginday": None,
                            "vac1-phase2-endday": None,
                            "vac1-frac-10": None,
                            "vac1-frac-20": None,
                            "vac1-frac-30": None,
                            "vac1-frac-40": None,
                            "vac1-frac-50": None,
                            "vac1-frac-60": None,
                            "vac1-frac-70": None,
                            "vac1-frac-80": None,
                            "vac2-protect-duration": None,
                            "vac2-efficacy-halflife": None,
                            "vac2-efficacy-slope": None,
                            "vac2-phase1-beginday": None,
                            "vac2-phase1-endday": None,
                            "vac2-phase1-dpd": None,
                            "vac2-phase2-beginday": None,
                            "vac2-phase2-endday": None,
                            "vac2-frac-10": None,
                            "vac2-frac-20": None,
                            "vac2-frac-30": None,
                            "vac2-frac-40": None,
                            "vac2-frac-50": None,
                            "vac2-frac-60": None,
                            "vac2-frac-70": None,
                            "vac2-frac-80": None}
        ## base param_sim with suggested default values
        self.param_sim_default =   {"introday": -1,
                                    "tf": 365,
                                    "symp-frac": 0.25,
                                    "dev-hosp-young": 1.0,
                                    "dev-hosp-mid": 0.5,
                                    "dev-hosp-old": 0.5,
                                    "dev-icu-frac": 1.0,
                                    "dev-icu-frac-phase2": 1.0,
                                    "dev-icu-frac-phase2beginday": 150,
                                    "dev-len-hospstay": 1.0,
                                    "dev-ventdeath-mid": 0.7,
                                    "rel-beta-hosp": 0.2,
                                    "beta": [1.3],
                                    "death-prob-home-60": 0.01,
                                    "death-prob-home-70": 0.03,
                                    "death-prob-home-80": 0.05,
                                    "death-prob-home-80": 0.05,
                                    "earlymarch-hosp-factor": 1.0,
                                    "hosp-frac-10": 0.02,
                                    "hosp-frac-20": 0.02,
                                    "hosp-frac-30": 0.04,
                                    "hosp-frac-40": 0.06,
                                    "hosp-frac-50": 0.15,
                                    "hosp-frac-60": 0.20,
                                    "hosp-frac-70": 0.25,
                                    "hosp-frac-80": 0.35,
                                    "mean-time-vent": 10.8,
                                    "prob-icu-vent": 0.75,
                                    "self-isolation-factor": 1.0,
                                    "susc-0-20": 0.6,
                                    "susc-60-100": 1.0,
                                    "time-symp-to-hosp": 7.0,
                                    "min-mixinglevel-00": 0.0,
                                    "min-mixinglevel-10": 0.0,
                                    "min-mixinglevel-20": 0.0,
                                    "min-mixinglevel-30": 0.0,
                                    "min-mixinglevel-40": 0.0,
                                    "min-mixinglevel-50": 0.0,
                                    "min-mixinglevel-60": 0.0,
                                    "min-mixinglevel-70": 0.0,
                                    "min-mixinglevel-80": 0.0,
                                    "contact-rate-10": 1.0,
                                    "contact-rate-20": 1.0,
                                    "contact-rate-30": 1.0,
                                    "contact-rate-40": 1.0,
                                    "contact-rate-50": 1.0,
                                    "contact-rate-60": 1.0,
                                    "contact-rate-70": 1.0,
                                    "contact-rate-80": 1.0,
                                    "firstlockdown-endday": 200,
                                    "contact-rate-postld-10": 1.0,
                                    "contact-rate-postld-20": 1.0,
                                    "contact-rate-postld-30": 1.0,
                                    "contact-rate-postld-40": 1.0,
                                    "contact-rate-postld-50": 1.0,
                                    "contact-rate-postld-60": 1.0,
                                    "contact-rate-postld-70": 1.0,
                                    "contact-rate-postld-80": 1.0,
                                    "contact-coeff-00": 1.0,
                                    "contact-coeff-10": 1.0,
                                    "contact-coeff-20": 1.0,
                                    "contact-coeff-30": 1.0,
                                    "contact-coeff-40": 1.0,
                                    "contact-coeff-50": 1.0,
                                    "contact-coeff-60": 1.0,
                                    "contact-coeff-70": 1.0,
                                    "contact-coeff-80": 1.0,
                                    "contact-coeff-postld-00": 1.0,
                                    "contact-coeff-postld-10": 1.0,
                                    "contact-coeff-postld-20": 1.0,
                                    "contact-coeff-postld-30": 1.0,
                                    "contact-coeff-postld-40": 1.0,
                                    "contact-coeff-postld-50": 1.0,
                                    "contact-coeff-postld-60": 1.0,
                                    "contact-coeff-postld-70": 1.0,
                                    "contact-coeff-postld-80": 1.0,
                                    "normalcy-beginday": 367,
                                    "contact-coeff-norm-00": 1.0,
                                    "contact-coeff-norm-10": 1.0,
                                    "contact-coeff-norm-20": 1.0,
                                    "contact-coeff-norm-30": 1.0,
                                    "contact-coeff-norm-40": 1.0,
                                    "contact-coeff-norm-50": 1.0,
                                    "contact-coeff-norm-60": 1.0,
                                    "contact-coeff-norm-70": 1.0,
                                    "contact-coeff-norm-80": 1.0,
                                    "len-nat-immunity": 540,
                                    "vac1-protect-duration": 0.1,
                                    "vac1-efficacy-halflife": 0.1,
                                    "vac1-efficacy-slope": 50.0,
                                    "vac1-phase1-beginday": 366,
                                    "vac1-phase1-endday": 367,
                                    "vac1-phase1-dpd": 0,
                                    "vac1-phase2-beginday": 368,
                                    "vac1-phase2-endday": 369,
                                    "vac1-frac-10": 0.0,
                                    "vac1-frac-20": 0.0,
                                    "vac1-frac-30": 0.0,
                                    "vac1-frac-40": 0.0,
                                    "vac1-frac-50": 0.0,
                                    "vac1-frac-60": 0.0,
                                    "vac1-frac-70": 0.0,
                                    "vac1-frac-80": 0.0,
                                    "vac2-protect-duration": 0.1,
                                    "vac2-efficacy-halflife": 0.1,
                                    "vac2-efficacy-slope": 50.0,
                                    "vac2-phase1-beginday": 340,
                                    "vac2-phase1-endday": 341,
                                    "vac2-phase1-dpd": 0,
                                    "vac2-phase2-beginday": 342,
                                    "vac2-phase2-endday": 343,
                                    "vac2-frac-10": 0.0,
                                    "vac2-frac-20": 0.0,
                                    "vac2-frac-30": 0.0,
                                    "vac2-frac-40": 0.0,
                                    "vac2-frac-50": 0.0,
                                    "vac2-frac-60": 0.0,
                                    "vac2-frac-70": 0.0,
                                    "vac2-frac-80": 0.0}
        ## tuple showing the valid value range (min,max) of params
        ## e.g. (0,1) means the param value must be in [0,1]
        ## to be used in MCMC log_prior calculation
        ## last updated: 2020-06-01 (as of odesim version 5)
        self.param_sim_value_range =   {"introday": (None,None),
                                        "tf": (1,None),
                                        "symp-frac": (0.1, 0.3),
                                        "dev-hosp-young": (0.5, 2.5),
                                        "dev-hosp-mid": (0.1, 2.0),
                                        "dev-hosp-old": (0.1, 1.2),
                                        "dev-icu-frac": (0.5, 1.5),
                                        "dev-icu-frac-phase2": (0.5, 1.5),
                                        "dev-icu-frac-phase2beginday": (),
                                        "dev-len-hospstay": (0.5, 2.0),
                                        "dev-ventdeath-mid": (0.5, 1.3),
                                        "rel-beta-hosp": (0.0, 0.3),
                                        "beta": (0.0, None),
                                        "death-prob-home-60": (0.0, 0.1),
                                        "death-prob-home-70": (0.0, 0.2),
                                        "death-prob-home-80": (0.0, 0.3),
                                        "death-prob-home-80": (0.02, 0.3),
                                        "earlymarch-hosp-factor": (0.0, 1.0),
                                        "hosp-frac-10": (0.0, 0.2),
                                        "hosp-frac-20": (0.0, 0.2),
                                        "hosp-frac-30": (0.0, 0.2),
                                        "hosp-frac-40": (0.0, 0.5),
                                        "hosp-frac-50": (0.0, 0.5),
                                        "hosp-frac-60": (0.0, 1.0),
                                        "hosp-frac-70": (0.0, 1.0),
                                        "hosp-frac-80": (0.0, 1.0),
                                        "mean-time-vent": (0.0, None),
                                        "prob-icu-vent": (0.5, 0.99999),
                                        "self-isolation-factor": (0.0, 1.0),
                                        "susc-0-20": (0.2, 1.5),
                                        "susc-60-100": (0.5, 1.5),
                                        "time-symp-to-hosp": (3.0, 10.0),
                                        "min-mixinglevel-00": (0.0, 1.0),
                                        "min-mixinglevel-10": (0.0, 1.0),
                                        "min-mixinglevel-20": (0.0, 1.0),
                                        "min-mixinglevel-30": (0.0, 1.0),
                                        "min-mixinglevel-40": (0.0, 1.0),
                                        "min-mixinglevel-50": (0.0, 1.0),
                                        "min-mixinglevel-60": (0.0, 1.0),
                                        "min-mixinglevel-70": (0.0, 1.0),
                                        "min-mixinglevel-80": (0.0, 1.0),
                                        "contact-rate-10": (0.2 ,4.0),
                                        "contact-rate-20": (0.2 ,4.0),
                                        "contact-rate-30": (0.2 ,4.0),
                                        "contact-rate-40": (0.2 ,4.0),
                                        "contact-rate-50": (0.2 ,4.0),
                                        "contact-rate-60": (0.2 ,4.0),
                                        "contact-rate-70": (0.2 ,4.0),
                                        "contact-rate-80": (0.2 ,4.0),
                                        "firstlockdown-endday": (100, 250),
                                        "contact-rate-postld-10": (0.2, 4.0),
                                        "contact-rate-postld-20": (0.2, 4.0),
                                        "contact-rate-postld-30": (0.2, 4.0),
                                        "contact-rate-postld-40": (0.2, 4.0),
                                        "contact-rate-postld-50": (0.2, 4.0),
                                        "contact-rate-postld-60": (0.2, 4.0),
                                        "contact-rate-postld-70": (0.2, 4.0),
                                        "contact-rate-postld-80": (0.2, 4.0),
                                        "contact-coeff-00": (0.0, None),
                                        "contact-coeff-10": (0.0, None),
                                        "contact-coeff-20": (0.0, None),
                                        "contact-coeff-30": (0.0, None),
                                        "contact-coeff-40": (0.0, None),
                                        "contact-coeff-50": (0.0, None),
                                        "contact-coeff-60": (0.0, None),
                                        "contact-coeff-70": (0.0, None),
                                        "contact-coeff-80": (0.0, None),
                                        "contact-coeff-postld-00": (0.0, None),
                                        "contact-coeff-postld-10": (0.0, None),
                                        "contact-coeff-postld-20": (0.0, None),
                                        "contact-coeff-postld-30": (0.0, None),
                                        "contact-coeff-postld-40": (0.0, None),
                                        "contact-coeff-postld-50": (0.0, None),
                                        "contact-coeff-postld-60": (0.0, None),
                                        "contact-coeff-postld-70": (0.0, None),
                                        "contact-coeff-postld-80": (0.0, None),
                                        "normalcy-beginday": (61, None),
                                        "contact-coeff-norm-00": (0.0, None),
                                        "contact-coeff-norm-10": (0.0, None),
                                        "contact-coeff-norm-20": (0.0, None),
                                        "contact-coeff-norm-30": (0.0, None),
                                        "contact-coeff-norm-40": (0.0, None),
                                        "contact-coeff-norm-50": (0.0, None),
                                        "contact-coeff-norm-60": (0.0, None),
                                        "contact-coeff-norm-70": (0.0, None),
                                        "contact-coeff-norm-80": (0.0, None),
                                        "len-nat-immunity": (0.0, None),
                                        "vac1-protect-duration": (0.1, None),
                                        "vac1-efficacy-halflife": (0.1, None),
                                        "vac1-efficacy-slope": (0.0, None),
                                        "vac1-phase1-beginday": (61, None),
                                        "vac1-phase1-endday": (62, None),
                                        "vac1-phase1-dpd": (0, None),
                                        "vac1-phase2-beginday": (63, None),
                                        "vac1-phase2-endday": (64, None),
                                        "vac1-frac-10": (0.0, 1.0),
                                        "vac1-frac-20": (0.0, 1.0),
                                        "vac1-frac-30": (0.0, 1.0),
                                        "vac1-frac-40": (0.0, 1.0),
                                        "vac1-frac-50": (0.0, 1.0),
                                        "vac1-frac-60": (0.0, 1.0),
                                        "vac1-frac-70": (0.0, 1.0),
                                        "vac1-frac-80": (0.0, 1.0),
                                        "vac2-protect-duration": (0.1, None),
                                        "vac2-efficacy-halflife": (0.1, None),
                                        "vac2-efficacy-slope": (0.0, None),
                                        "vac2-phase1-beginday": (65, None),
                                        "vac2-phase1-endday": (66, None),
                                        "vac2-phase1-dpd": (0, None),
                                        "vac2-phase2-beginday": (67, None),
                                        "vac2-phase2-endday": (68,None),
                                        "vac2-frac-10": (0.0, 1.0),
                                        "vac2-frac-20": (0.0, 1.0),
                                        "vac2-frac-30": (0.0, 1.0),
                                        "vac2-frac-40": (0.0, 1.0),
                                        "vac2-frac-50": (0.0, 1.0),
                                        "vac2-frac-60": (0.0, 1.0),
                                        "vac2-frac-70": (0.0, 1.0),
                                        "vac2-frac-80": (0.0, 1.0)
                                        }

        self.current_traj = None
        self.indices = self.Get_Indices()
        self.version = self.Get_Version()
    
    ## get odesim version
    ## to be used in setting command
    def Get_Version(self):
        outs = run([str(self.cpp_dir) + "/odesim", "none", "-version"], capture_output=True).stdout
        version = outs.decode("utf-8").split('\n')[0]
        return int(version)
        

    ## get indices of compartments from the simulation
    ## to be used in getting compartment values
    def Get_Indices(self):
        outs = run([str(self.cpp_dir) + "/odesim", "none", "-printIndices"], capture_output=True).stdout
        ## dict( zip( [key_list], [value_list]) )
        indices = dict(
            zip([
                x.split(" ")[0] for x in outs.decode("utf-8").split('\n')[:-1]
            ], [
                int(x.split(" ")[1])
                for x in outs.decode("utf-8").split('\n')[:-1]
            ]))
        return indices

    ## construct command to run the simulation
    def Set_Sim_CMD_str(self,
                    output_file=None,
                    loc=None,
                    contact_matrix=False,
                    steps_per_day=None,
                    binary_output=False,
                    death_rate_only=False,
                    symp_frac_equal=None,
                    symp_frac_davies=False,
                    **kwargs):
        '''
        output_file: str
            save odesim output to output_file (instead of standard output)
        loc: str
            either "RI", "MA", or "PA"
        contact_matrix: bool
            turning on using social contact matrix, available from odesim v5
        steps_per_day: int
            the number of integration steps per day, available from version 5
        binary_output: bool
            turning on printing binary output (instead of text)
        kwargs: dict
            command-line parameters for odesim
            keywords must match those in param_dict
        '''

        cmd = str(self.cpp_dir) + "/odesim"  # base command

        ## construct full command with arguments
        if output_file is not None: cmd = " ".join([cmd, output_file])
        else: cmd = " ".join([cmd, "none"])
        ## contact matrix is available from (?) odesim version 5
        if contact_matrix and self.version > 4: cmd = " ".join([cmd, "-scm"])
        ## location: either "RI", "MA", or "PA"
        if loc is not None:
            if loc not in ["RI", "MA", "PA"] : raise ValueError('Location must be either "RI", "MA" or "PA".')
            else: cmd = " ".join([cmd, "-loc", loc])
        ## integration steps per day
        if steps_per_day is not None and self.version > 4:
            if int(steps_per_day) < 2: raise ValueError('steps_per_day argument should be higher than 1.')
            else: cmd = " ".join([cmd, "-steps-per-day", str(steps_per_day)])
        ## tell odesim to produce binary output, available from version 5
        if binary_output and self.version > 4:
            cmd = " ".join([cmd, "-binary-output"])
        ## tell odesim to produce death-rate-output only, available from version 5
        if death_rate_only and self.version > 4:
            cmd = " ".join([cmd, "-death-rate-output"])
        ## tell odesim to use symp-frac-equal, available from version 5
        if symp_frac_equal is not None and self.version > 4:
            cmd = " ".join([cmd, "-symp-frac-equal", str(symp_frac_equal)])
        ## tell odesim to use symp-frac-davies, available from version 5
        if symp_frac_davies and self.version > 4:
            cmd = " ".join([cmd, "-symp-frac-davies"])   

        ## handle **kwargs
        beta_str = ""
        for k,v in kwargs.items():
            if k in self.param_sim_version:
                valid_version = self.param_sim_version[k]
                if self.version >= valid_version[0] and self.version < valid_version[1] and v is not None:
                    if k == "beta": beta_str = " ".join(["-beta"] + [str(x) for x in v])
                    else: cmd = " ".join([cmd, "-"+k, str(v)])
        cmd = " ".join([cmd, beta_str])

        return cmd
    
    ## construct command to run the simulation
    def Set_Sim_CMD(self,
                    output_file=None,
                    loc=None,
                    contact_matrix=False,
                    steps_per_day=None,
                    binary_output=False,
                    death_rate_only=False,
                    symp_frac_equal=None,
                    symp_frac_davies=False,
                    **kwargs):
        '''
        output_file: str
            save odesim output to output_file (instead of standard output)
        loc: str
            either "RI", "MA", or "PA"
        contact_matrix: bool
            turning on using social contact matrix, available from odesim v5
        steps_per_day: int
            the number of integration steps per day, available from version 5
        binary_output: bool
            turning on printing binary output (instead of text)
        kwargs: dict
            command-line parameters for odesim
            keywords must match those in param_dict
        '''

        cmd = [str(self.cpp_dir) + "/odesim"]  # base command

        ## construct full command with arguments
        if output_file is not None: cmd.append(output_file)
        else: cmd.append("none")
        ## contact matrix is available from odesim version 5
        if contact_matrix and self.version > 4: cmd.append("-scm")
        ## location: either "RI", "MA", or "PA"
        if loc is not None:
            if loc not in ["RI", "MA", "PA"] : raise ValueError('Location must be either "RI", "MA" or "PA".')
            else: cmd.extend(["-loc", loc])
        ## integration steps per day
        if steps_per_day is not None and self.version > 4:
            if int(steps_per_day) < 2.: raise ValueError('steps_per_day argument should be higher than 1.')
            else: cmd.extend(["-steps-per-day", str(steps_per_day)])
        ## tell odesim to produce binary output, available from version 5
        if binary_output and self.version > 4:
            cmd.append("-binary-output")
        ## tell odesim to produce death-rate-output only, available from version 5
        if death_rate_only and self.version > 4:
            cmd.append("-death-rate-output")
        ## tell odesim to use symp-frac-equal, available from version 5
        if symp_frac_equal is not None and self.version > 4:
            cmd.extend(["-symp-frac-equal", str(symp_frac_equal)])
        ## tell odesim to use symp-frac-davies, available from version 5
        if symp_frac_davies and self.version > 4:
            cmd.append("-symp-frac-davies")

        
        ## handle **kwargs
        beta_ls = []
        for k,v in kwargs.items():
            if k in self.param_sim_version:
                valid_version = self.param_sim_version[k]
                if self.version >= valid_version[0] and self.version < valid_version[1] and v is not None:
                    if k == "beta":
                        beta_ls.append("-beta")
                        beta_ls.extend([str(x) for x in v])
                    else:
                        cmd.extend(["-"+k, str(v)])
        cmd.extend(beta_ls)

        return cmd

    ## calling Popen to run the simulation
    def Run_Sim_CMD(self, output_file=None, loc=None, contact_matrix=False, 
                    steps_per_day=None, binary_output=False, death_rate_only=False,
                    symp_frac_equal=None, symp_frac_davies=False, print_cmd=False,
                    param_dict=dict()):
        ## unset contact_matrix if version < 5
        if contact_matrix and self.version < 5: contact_matrix = False
        ## unset steps_per_day if version < 5
        if steps_per_day is not None and self.version < 5: steps_per_day = None
        ## unset binary_output if version < 5
        if binary_output and self.version < 5: binary_output = False
        ## construct command
        # cmd = self.Set_Sim_CMD_str(output_file=output_file, loc=loc, contact_matrix=contact_matrix, 
        #                            steps_per_day=steps_per_day, binary_output=binary_output,
        #                            **param_dict)
        cmd = self.Set_Sim_CMD(output_file=output_file, loc=loc, contact_matrix=contact_matrix, 
                               steps_per_day=steps_per_day, binary_output=binary_output, death_rate_only=death_rate_only,
                               symp_frac_equal=symp_frac_equal, symp_frac_davies=symp_frac_davies,
                               **param_dict)
        if print_cmd:
            print(cmd)

        ## run command
        if output_file:
            proc = run(cmd, cwd=Path.cwd())
            if proc.returncode == 0:
                if binary_output:
                    if self.version < 6: self.current_traj = fromfile(Path(Path.cwd(), output_file)).reshape((-1, self.indices['STARTK'] + self.indices['NUMAC'] + 1))
                    else: self.current_traj = fromfile(Path(Path.cwd(), output_file)).reshape((-1, self.indices['DIMENSION'] + 1 ))
                else: self.current_traj = genfromtxt(Path(Path.cwd(), output_file), delimiter="\t")
            else:
                raise Exception("odesim failed to run")
        else:
            proc = run(cmd, capture_output=True)
            if proc.returncode == 0:
                if binary_output: 
                    if self.version < 6: self.current_traj = frombuffer(proc.stdout).reshape((-1,self.indices['STARTK'] + self.indices['NUMAC'] + 1))
                    else: self.current_traj = frombuffer(proc.stdout).reshape((-1,self.indices['DIMENSION'] + 1))
                else: self.current_traj = genfromtxt(BytesIO(proc.stdout), delimiter="\t")
            else: 
                raise Exception("odesim failed to run")
        
        return self.current_traj


    def Get_Compartment_START_Index(self, compartment):
        '''
        compartment: str
            odesim v4 accepts "E", "I", "A", "HA", "CA", "V", "CR", "HR", "D", "R", "J", "K"
            odesim v5 accepts all from v4 and "DHOSP", "RHOSP"
        '''
        key = "START{}".format(compartment)
        if key in self.indices:
            return self.indices[key]
        else:
            raise Exception('START{} is not available in this odesim version.'.format(compartment))
            # )

    def Get_Compartment_NUM_Stages(self, compartment):
        '''
        compartment: str
            odesim v4 accepts "E", "I", "A", "HA", "CA", "V", "CR", "HR", "D", "R", "J", "K"
            odesim v5 accepts all from v4 and "DHOSP", "RHOSP"
        '''
        key = "START{}".format(compartment)
        key_num = "NUM{}".format(compartment)
        if key in self.indices:
            if key_num in self.indices: return self.indices[key_num]
            else: return 1
        else:
            raise Exception('NUM{} is not available in this odesim version.'.format(compartment))


    def Get_Compartment_Age_Group(self, compartment, age_group_index, traj=None):
        '''
        traj: (ndays, ncolumns) numpy array where
            ndays: total number of simulated days
            ncolumns: the number of columns printed by odesim simulation at each timestep
        compartment: str
            odesim v4 accepts "E", "I", "A", "HA", "CA", "V", "CR", "HR", "D", "R", "J", "K"
            odesim v5 accepts all from v4 and "DHOSP", "RHOSP"
        age_group_index: int 
            range from 0 to (odesim's "NUMAC" - 1)
        '''
        if traj is None and self.current_traj is None:
            raise Exception("Trajectory is missing.")
        if age_group_index < 0 or age_group_index >= self.indices["NUMAC"]:
            raise ValueError("age_group_index must be in range 0-{}".format(self.indices["NUMAC"] - 1))
        if traj is None:
            traj = self.current_traj
        offset = 1  ## the first column in sim output is "day"
        comp_i = self.Get_Compartment_START_Index(compartment)
        num_i = self.Get_Compartment_NUM_Stages(compartment)
        idx = [offset + comp_i + self.indices["NUMAC"]*stage_i + age_group_index for stage_i in range(num_i)]

        # return traj[:, offset + comp_i + age_group_index]
        if num_i == 1: return traj[:, idx[0]]
        else: return npsum(traj[:, idx], axis=1)

    def Get_Compartment_All_Ages(self, compartment, traj=None):
        '''
        traj: (ndays, ncolumns) numpy array where
            ndays: total number of simulated days
            ncolumns: the number of columns printed by odesim simulation at each timestep
        compartment: str
            odesim v4 accepts "E", "I", "A", "HA", "CA", "V", "CR", "HR", "D", "R", "J", "K"
            odesim v5 accepts all from v4 and "DHOSP", "RHOSP"
        '''
        if traj is None and self.current_traj is None:
            raise Exception("Trajectory is missing.")
        if traj is None:
            traj = self.current_traj
        offset = 1  ## the first column in sim output is "day"
        comp_i = self.Get_Compartment_START_Index(compartment)
        num_i = self.Get_Compartment_NUM_Stages(compartment)

        return npsum(traj[:, offset + comp_i : offset + comp_i + self.indices["NUMAC"]*num_i ],
                     axis=1)

# %%
################
if __name__ == "__main__":

    # %%
    sim = COVID_SIM(Path("~/Code/covid19/cpp-v5-discharges-nonhospdeaths/")) ## cpp-v4-6e-severe-classes/")) ## 
    print("version", sim.version)
    print(sim.indices)

    # %%
    sim.param_sim["tf"] = 100
    sim.param_sim["beta"] = [2.3, 1.9, 1.4]
    sim.param_sim["hosp-frac-70"] = 0.9 ## version 5
    sim.Run_Sim_CMD(output_file=None, loc=None, contact_matrix=False, 
                    steps_per_day=5, binary_output=True,
                    param_dict=sim.param_sim)
    # print(sim.Get_Compartment_NUM_Stages(compartment="I"))
    print(sim.current_traj.shape)
    print(sim.Get_Compartment_Age_Group(compartment="V", age_group_index=7))
    print(sim.Get_Compartment_All_Ages(compartment="RHOSP")) ## version 5 


# %%
