#!/bin/bash

## usage:
## ./run_tf1800.sh 
##                  -r odesim_directory;                default: current directory
##                  -o output_directory;                default: output
##                  -s campaign_start_day, can be a space-delimited list;       default: 426
##                  -e campaign_end_day, can be a space-delimited list;         default: 467
##                  -d dose_per_day, can be a space-delimited list;             default: 10000 
##                  -i natural_immunity;                default: 540
##                  -c normalcy_startday;               default: 1098
##                  -u vaccine_protection_duration;     default: 540
##                  -t vaccine_halflife;                default: 180
##                  -p vaccine_slope;                   default: 2
## NOTE:
##      1. campaign_start_day, campaign_end_day, dose_per_day, if given as a list, must have equal number of elements
##      2. -tf option for odesim = 60 + the number of betas given in the "beta" string below

odesim_dir='.'
outdir='output'
startday=426
endday=467
dpd=10000
nat_imm=540
normalcyday=1098
vac_duration=540
vac_halflife=180
vac_slope=2

## other odesim params 

default_config="-no-test-before-vaccinate -loc RI -introday 55 -steps-per-day 2 -symp-frac-davies -contact-coeff-00 3.3738859 -contact-coeff-10 3.868847 -contact-coeff-20 4.6584151 -contact-coeff-30 5.1108777 -contact-coeff-40 4.2338707 -contact-coeff-50 4.0706229 -contact-coeff-60 2.5729902 -contact-coeff-70 3.7762 -contact-coeff-80 12.5471046 -contact-coeff-postld-00 1.0934778 -contact-coeff-postld-10 2.1349499 -contact-coeff-postld-20 0.729322 -contact-coeff-postld-30 0.7973113 -contact-coeff-postld-40 0.5223126 -contact-coeff-postld-50 0.4670234 -contact-coeff-postld-60 0.3283015 -contact-coeff-postld-70 0.401475 -contact-coeff-postld-80 1.320971 -firstlockdown-endday 142.2869914 -mean-time-vent 7.9725535 -death-prob-home-60 0.0163857 -death-prob-home-70 0.069957 -death-prob-home-80 0.2256238 -time-symp-to-hosp 2.1514175 -dev-len-hospstay 0.4875795 -dev-icu-frac 0.4047451 -dev-icu-frac-phase2 0.3543851 -dev-icu-frac-phase2beginday 157.3115232 -prob-icu-vent 0.6754456 -dev-ventdeath-mid 0.8501806 -hosp-frac-10 0.0176092 -hosp-frac-20 0.0251054 -hosp-frac-30 0.0362311 -hosp-frac-40 0.0602989 -hosp-frac-50 0.0873153 -hosp-frac-60 0.1629561 -hosp-frac-70 0.2926905 -hosp-frac-80 0.2634612"


# trough/peak = 0.1
beta="1.6727976 1.3091595 1.0556721 0.839467 0.6612355 0.5318107 0.4400402 0.3742544 0.3336967 0.2996406 0.2729291 0.2547003 0.2390196 0.2287049 0.2211882 0.2135892 0.2062044 0.1985657 0.1948971 0.1892483 0.1783025 0.1718373 0.1674553 0.1618883 0.1591671 0.1565791 0.155222 0.1547845 0.1531179 0.1516463 0.1469013 0.1414765 0.1356835 0.1281881 0.1214586 0.1142545 0.1079651 0.0998169 0.0942869 0.0902883 0.0875188 0.0861613 0.0851119 0.0838469 0.0846169 0.0855279 0.0855725 0.0856781 0.0851056 0.0831568 0.0803501 0.0763618 0.0727229 0.0699567 0.066197 0.0623294 0.0591034 0.0552461 0.0524101 0.048635 0.0465685 0.0458828 0.0468066 0.0480444 0.0503944 0.0528407 0.0558019 0.0569346 0.0581211 0.0588087 0.0593474 0.058193 0.0566061 0.056194 0.0571379 0.0590063 0.0616593 0.0634172 0.0676003 0.07469 0.0857409 0.0965157 0.113263 0.1344348 0.1585997 0.1863518 0.2260575 0.2611772 0.2945457 0.3283191 0.357635 0.3868245 0.4093747 0.4175377 0.4246134 0.4265253 0.4129889 0.4171379 0.4038737 0.3918221 0.3705054 0.3500103 0.3305635 0.3021442 0.293615 0.2845611 0.2759113 0.2699475 0.2649497 0.2712738 0.2789295 0.2855462 0.3019285 0.317321 0.3335838 0.3551688 0.3783127 0.3904646 0.4084734 0.4188742 0.4410836 0.4614731 0.4687238 0.4800018 0.4939878 0.487615 0.5064311 0.5060004 0.5104718 0.5243046 0.534116 0.5293655 0.5136221 0.5030991 0.4887424 0.4907996 0.4883806 0.4959582 0.4964565 0.5059589 0.4983837 0.49302 0.4964822 0.4996463 0.4999052 0.4953054 0.5011131 0.5030477 0.5004874 0.4996323 0.4951594 0.4832767 0.4682172 0.4563955 0.4450037 0.4228513 0.4111123 0.395628 0.3893881 0.3847026 0.3812418 0.3753286 0.3695844 0.3793789 0.383356 0.3892003 0.385007 0.386506 0.3716624 0.3765839 0.3759507 0.3574669 0.3559739 0.3577437 0.3603862 0.3612102 0.3675355 0.3735665 0.3779986 0.3865911 0.3802716 0.3724093 0.3719672 0.3857572 0.4022878 0.4190901 0.4388562 0.4524915 0.4739582 0.4964955 0.5164054 0.5358147 0.5529294 0.5704789 0.5778226 0.5936818 0.6159678 0.6169427 0.6127916 0.6235248 0.617918 0.6025216 0.590118 0.5966963 0.5890876 0.584956 0.5800465 0.600875 0.6192757 0.6417883 0.6396446 0.6353462 0.6359176 0.6346613 0.6360578 0.6331754 0.6370787 0.6340793 0.6227533 0.6285719 0.6225454 0.6215102 0.6154268 0.6124175 0.6120027 0.6100525 0.6101937 0.6151179 0.6284947 0.6270344 0.6247281 0.623903 0.6204963 0.6118679 0.6098525 0.6087834 0.6024977 0.6051303 0.6081419 0.6019591 0.6034065 0.5993969 0.6004755 0.599124 0.6028464 0.6044168 0.6131865 0.6183919 0.6235973 0.6235772 0.6256727 0.6268633 0.626969 0.631818 0.6281357 0.6239917 0.6155393 0.6038295 0.5981786 0.5865631 0.5727228 0.5607739 0.5492601 0.5423473 0.5444192 0.54388 0.5451427 0.545438 0.5668876 0.5975211 0.6160267 0.6424269 0.6639603 0.6920593 0.7132619 0.7379502 0.7846868 0.7791358 0.8053503 0.8349254 0.8274932 0.8274932 0.8133633 0.7992334 0.7851035 0.7709736 0.7568437 0.7427138 0.7285839 0.7144539 0.700324 0.6861941 0.6720642 0.6579343 0.6438044 0.6296745 0.6155446 0.6014147 0.5872848 0.5731549 0.559025 0.5448951 0.5307652 0.5166353 0.5025054 0.4883755 0.4742456 0.4601157 0.4459857 0.4318558 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.4177259 0.418445 0.4191638 0.4198822 0.4206 0.4213169 0.4220328 0.4227473 0.4234604 0.4241718 0.4248813 0.4255887 0.4262937 0.4269962 0.4276959 0.4283927 0.4290863 0.4297766 0.4304633 0.4311462 0.4318251 0.4324999 0.4331703 0.4338361 0.4344971 0.4351532 0.4358041 0.4364496 0.4370896 0.4377239 0.4383522 0.4389744 0.4395904 0.4401998 0.4408026 0.4413985 0.4419875 0.4425692 0.4431436 0.4437104 0.4442696 0.4448209 0.4453641 0.4458992 0.4464259 0.4469441 0.4474537 0.4479544 0.4484462 0.4489289 0.4494023 0.4498664 0.4503209 0.4507658 0.4512009 0.451626 0.4520411 0.4524461 0.4528408 0.453225 0.4535987 0.4539618 0.4543142 0.4546557 0.4549863 0.4553059 0.4556143 0.4559114 0.4561973 0.4564718 0.4567347 0.4569862 0.457226 0.457454 0.4576703 0.4578748 0.4580674 0.458248 0.4584166 0.4585732 0.4587177 0.45885 0.4589701 0.459078 0.4591736 0.459257 0.4593281 0.4593868 0.4594332 0.4594672 0.4594889 0.4594981 0.4594951 0.4594796 0.4594517 0.4594115 0.459359 0.4592941 0.4592169 0.4591273 0.4590256 0.4589116 0.4587853 0.4586469 0.4584964 0.4583338 0.4581592 0.4579726 0.4577741 0.4575637 0.4573415 0.4571075 0.4568619 0.4566047 0.456336 0.4560558 0.4557643 0.4554615 0.4551475 0.4548224 0.4544863 0.4541394 0.4537816 0.4534132 0.4530342 0.4526447 0.4522449 0.4518349 0.4514147 0.4509846 0.4505446 0.4500949 0.4496355 0.4491668 0.4486887 0.4482015 0.4477052 0.4472"


## vaccination stratgies:
##      0. no vaccine
##      1. random: "-vac1-ratio-10 0.11111111 -vac1-ratio-20 0.11111111 -vac1-ratio-30 0.11111111 -vac1-ratio-40 0.11111111 -vac1-ratio-50 0.11111111 -vac1-ratio-60 0.11111111 -vac1-ratio-70 0.11111111 -vac1-ratio-80 0.11111111" 
##      2. 16-29 only*: "-vac1-ratio-10 0.28571400 -vac1-ratio-20 0.71428500 -vac1-ratio-30 0.00000014 -vac1-ratio-40 0.00000014 -vac1-ratio-50 0.00000014 -vac1-ratio-60 0.00000014 -vac1-ratio-70 0.00000014 -vac1-ratio-80 0.00000014"
##      3. 30-59 only*: "-vac1-ratio-10 0.00000017 -vac1-ratio-20 0.00000017 -vac1-ratio-30 0.33333300 -vac1-ratio-40 0.33333300 -vac1-ratio-50 0.33333300 -vac1-ratio-60 0.00000017 -vac1-ratio-70 0.00000017 -vac1-ratio-80 0.00000017"
##      4. 60+ only*: "-vac1-ratio-10 0.00000017 -vac1-ratio-20 0.00000017 -vac1-ratio-30 0.00000017 -vac1-ratio-40 0.00000017 -vac1-ratio-50 0.00000017 -vac1-ratio-60 0.33333300 -vac1-ratio-70 0.33333300 -vac1-ratio-80 0.33333300"
##      5. random 50/50 excl. 30-59*: "-vac1-ratio-10 0.16666650 -vac1-ratio-20 0.16666650 -vac1-ratio-30 0.00000033 -vac1-ratio-40 0.00000033 -vac1-ratio-50 0.00000033 -vac1-ratio-60 0.16666650 -vac1-ratio-70 0.16666650 -vac1-ratio-80 0.16666650"
##      6. random 25/75 excl. 30-59*: "-vac1-ratio-10 0.08333325 -vac1-ratio-20 0.08333325 -vac1-ratio-30 0.00000033 -vac1-ratio-40 0.00000033 -vac1-ratio-50 0.00000033 -vac1-ratio-60 0.24999975 -vac1-ratio-70 0.24999975 -vac1-ratio-80 0.24999975"
##      7. random 75/25 excl. 30-59*: "-vac1-ratio-10 0.24999975 -vac1-ratio-20 0.24999975 -vac1-ratio-30 0.00000033 -vac1-ratio-40 0.00000033 -vac1-ratio-50 0.00000033 -vac1-ratio-60 0.08333325 -vac1-ratio-70 0.08333325 -vac1-ratio-80 0.08333325"
##      8. random among 20-49 and 70+*: "-vac1-ratio-10 0.00000025 -vac1-ratio-20 0.1999998 -vac1-ratio-30 0.1999998 -vac1-ratio-40 0.1999998 -vac1-ratio-50 0.00000025 -vac1-ratio-60 0.00000025 -vac1-ratio-70 0.1999998 -vac1-ratio-80 0.1999998"
##      9. random 50/50 among 20-49 and 70+*: "-vac1-ratio-10 0.00000025 -vac1-ratio-20 0.16666650 -vac1-ratio-30 0.16666650 -vac1-ratio-40 0.16666650 -vac1-ratio-50 0.00000025 -vac1-ratio-60 0.00000025 -vac1-ratio-70 0.24999975 -vac1-ratio-80 0.24999975"
##     10. random 25/75 among 20-49 and 70+*: "-vac1-ratio-10 0.00000025 -vac1-ratio-20 0.08333325 -vac1-ratio-30 0.08333325 -vac1-ratio-40 0.08333325 -vac1-ratio-50 0.00000025 -vac1-ratio-60 0.00000025 -vac1-ratio-70 0.374999625 -vac1-ratio-80 0.374999625"
##     11. random 75/25 among 20-49 and 70+*: "-vac1-ratio-10 0.00000025 -vac1-ratio-20 0.24999975 -vac1-ratio-30 0.24999975 -vac1-ratio-40 0.24999975 -vac1-ratio-50 0.00000025 -vac1-ratio-60 0.00000025 -vac1-ratio-70 0.124999875 -vac1-ratio-80 0.124999875"
declare -a sim=("-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 

"-vac1-ratio-10 0.2857142829 -vac1-ratio-20 0.7142857071 -vac1-ratio-30 0.0000000017 -vac1-ratio-40 0.0000000017 -vac1-ratio-50 0.0000000017 -vac1-ratio-60 0.0000000017 -vac1-ratio-70 0.0000000017 -vac1-ratio-80 0.0000000017" 
"-vac1-ratio-10 0.0000000008 -vac1-ratio-20 0.0000000023 -vac1-ratio-30 0.3333333300 -vac1-ratio-40 0.3333333300 -vac1-ratio-50 0.3333333300 -vac1-ratio-60 0.0000000023 -vac1-ratio-70 0.0000000023 -vac1-ratio-80 0.0000000023" 
"-vac1-ratio-10 0.0000000008 -vac1-ratio-20 0.0000000023 -vac1-ratio-30 0.0000000023 -vac1-ratio-40 0.0000000023 -vac1-ratio-50 0.0000000023 -vac1-ratio-60 0.3333333300 -vac1-ratio-70 0.3333333300 -vac1-ratio-80 0.3333333300"

"-vac1-frac-10 0 -vac1-frac-20 0.25 -vac1-frac-30 0.25 -vac1-frac-40 0 -vac1-frac-50 0 -vac1-frac-60 0.1666666667 -vac1-frac-70 0.1666666667 -vac1-frac-80 0.1666666667  
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 
"-vac1-frac-10 0 -vac1-frac-20 0.125 -vac1-frac-30 0.125 -vac1-frac-40 0 -vac1-frac-50 0 -vac1-frac-60 0.25 -vac1-frac-70 0.25 -vac1-frac-80 0.25  
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 
"-vac1-frac-10 0 -vac1-frac-20 0.375 -vac1-frac-30 0.375 -vac1-frac-40 0 -vac1-frac-50 0 -vac1-frac-60 0.0833333333 -vac1-frac-70 0.0833333333 -vac1-frac-80 0.0833333333  
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 

"-vac1-frac-10 0 -vac1-frac-20 0.2 -vac1-frac-30 0.2 -vac1-frac-40 0.2 -vac1-frac-50 0 -vac1-frac-60 0 -vac1-frac-70 0.2 -vac1-frac-80 0.2
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 

"-vac1-frac-10 0 -vac1-frac-20 0.16666667 -vac1-frac-30 0.16666667 -vac1-frac-40 0.16666667 -vac1-frac-50 0 -vac1-frac-60 0 -vac1-frac-70 0.25 -vac1-frac-80 0.25
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 

"-vac1-frac-10 0 -vac1-frac-20 0.08333333 -vac1-frac-30 0.08333333 -vac1-frac-40 0.08333333 -vac1-frac-50 0 -vac1-frac-60 0 -vac1-frac-70 0.375 -vac1-frac-80 0.375
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857" 

"-vac1-frac-10 0 -vac1-frac-20 0.25 -vac1-frac-30 0.25 -vac1-frac-40 0.25 -vac1-frac-50 0 -vac1-frac-60 0 -vac1-frac-70 0.125 -vac1-frac-80 0.125
-vac1-ratio-10 0.05 -vac1-ratio-20 0.1357142857 -vac1-ratio-30 0.1357142857 -vac1-ratio-40 0.1357142857 -vac1-ratio-50 0.1357142857 -vac1-ratio-60 0.1357142857 -vac1-ratio-70 0.1357142857 -vac1-ratio-80 0.1357142857")

## while there are still command line option to read in
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r|--simdir)
    odesim_dir="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outdir)
    outdir="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--startday)
    startday="$2"
    shift # past argument
    shift # past value
    ;;
    -e|--endday)
    endday="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--dpd)
    dpd="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--immunity)
    nat_imm="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--normalcyday)
    normalcyday="$2"
    shift # past argument
    shift # past value
    ;;
    -u|--vacduration)
    vac_duration="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--vachalflife)
    vac_halflife="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--vacslope)
    vac_slope="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    # POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

echo ""
echo "-r $odesim_dir -o $outdir -s $startday -e $endday -d $dpd -i $nat_imm -c $normalcyday -u $vac_duration -t $vac_halflife -p $vac_slope"




beta_a=($beta)
tf=$((60+${#beta_a[@]}))

echo "run 0 baseline no vaccine"
# ${odesim_dir}/./odesim ${outdir}/run-f${startday}_t${endday}_dpd${dpd}-0.tdl 
${odesim_dir}/./odesim ${outdir}/run_output_0.tdl -tf $tf $default_config -normalcy-beginday $normalcyday -len-nat-immunity $nat_imm -beta $beta


i=0 ## index of the last run before this loop
for vacfrac in "${sim[@]}"
do
    i=$(($i+1))
    echo "dpd $dpd run $i"
    # ${odesim_dir}/./odesim ${outdir}/run-f${startday:0:3}_t${endday:((${#endday}-3)):3}_dpd${dpd:0:5}-${i}.tdl 
    ${odesim_dir}/./odesim ${outdir}/run_output_${i}.tdl -tf $tf $default_config -normalcy-beginday $normalcyday -len-nat-immunity $nat_imm -vac1-protect-duration $vac_duration -vac1-efficacy-halflife $vac_halflife -vac1-efficacy-slope $vac_slope -vac1-phase1-beginday $startday -vac1-phase1-endday $endday -vac1-phase1-dpd $dpd $vacfrac -beta $beta
done
echo "total runs: $(($i+1))"

exit
