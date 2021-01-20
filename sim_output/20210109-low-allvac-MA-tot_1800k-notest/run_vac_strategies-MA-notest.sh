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
default_config="-binary-output -no-test-before-vaccinate -loc MA -symp-frac-davies -introday 55 -steps-per-day 2 -contact-coeff-00 3.2153261 -contact-coeff-10 0.6474441 -contact-coeff-20 2.9211862 -contact-coeff-30 3.4073948 -contact-coeff-40 3.2624102 -contact-coeff-50 2.9868068 -contact-coeff-60 2.1481057 -contact-coeff-70 3.2669315 -contact-coeff-80 13.7209337 -contact-coeff-postld-00 0.9358465 -contact-coeff-postld-10 9.5260095 -contact-coeff-postld-20 2.3152282 -contact-coeff-postld-30 2.7470483 -contact-coeff-postld-40 1.769564 -contact-coeff-postld-50 1.3999454 -contact-coeff-postld-60 1.0257234 -contact-coeff-postld-70 1.295245 -contact-coeff-postld-80 4.346264 -firstlockdown-endday 125.4204491 -mean-time-vent 5.4643665 -death-prob-home-70 0.0657149 -death-prob-home-80 0.3165232 -time-symp-to-hosp 9.9526423 -dev-len-hospstay 0.4896878 -dev-icu-frac 0.8760235 -dev-icu-frac-phase2 0.4053857 -dev-icu-frac-phase2beginday 157.3005705 -prob-icu-vent 0.7059609 -dev-ventdeath-mid 0.872809 -hosp-frac-10 0.0160262 -hosp-frac-20 0.0221684 -hosp-frac-30 0.0444339 -hosp-frac-40 0.0360453 -hosp-frac-50 0.099582 -hosp-frac-60 0.2361595 -hosp-frac-70 0.3845485 -hosp-frac-80 0.1526061"


# trough/peak = 0.1
beta="2.3746913 2.1638718 1.9832254 1.8092252 1.6767311 1.5211669 1.4298381 1.3059504 1.2046879 1.1032218 1.0172373 0.9258279 0.8230386 0.7445818 0.6616794 0.5880627 0.5290158 0.4748198 0.4225441 0.3749185 0.3364276 0.3031348 0.2800867 0.2614631 0.2466029 0.2367079 0.2316118 0.2265787 0.2255923 0.2188917 0.2088632 0.1983779 0.1816536 0.1662814 0.1473247 0.1276267 0.1110808 0.0939079 0.0804766 0.0713277 0.0664825 0.064596 0.0664336 0.0672511 0.0709947 0.0739457 0.0760269 0.0762482 0.0759755 0.076463 0.0751426 0.0722851 0.0703657 0.0679216 0.0666048 0.0654164 0.0652229 0.0649738 0.0655401 0.0659731 0.0664877 0.0665224 0.0659112 0.0645595 0.0641464 0.0627839 0.0598381 0.0559784 0.052044 0.046978 0.0417471 0.0364902 0.0317712 0.0274427 0.0242632 0.0222541 0.0218178 0.0228709 0.0258949 0.0300784 0.0357282 0.0422678 0.0490519 0.0555431 0.0616546 0.0678437 0.0735192 0.0776601 0.0819096 0.0842597 0.0841026 0.0830487 0.0795614 0.07512 0.0687972 0.0626376 0.0551261 0.0483213 0.0434319 0.040855 0.0399575 0.0408767 0.0425833 0.0454463 0.048695 0.052831 0.0585075 0.0633788 0.0705367 0.0751885 0.0797595 0.0844218 0.090564 0.0921595 0.0928834 0.0939008 0.0946791 0.0933223 0.0877528 0.0830776 0.0789428 0.0779193 0.0773143 0.0792376 0.0816392 0.0817419 0.0861357 0.0876823 0.0887811 0.0891417 0.0897708 0.0868865 0.0848263 0.0807699 0.0792984 0.0781015 0.0764933 0.0787482 0.0829532 0.0881596 0.0960992 0.1017504 0.1079787 0.1136491 0.1172905 0.1211117 0.1222918 0.1205102 0.11898 0.1154973 0.1110606 0.1075549 0.1040577 0.0999301 0.0965815 0.0933062 0.0904019 0.0885198 0.0877442 0.0870144 0.0849127 0.0842054 0.0823769 0.0805384 0.081433 0.082529 0.0821704 0.082529 0.0853797 0.0866034 0.0861529 0.0875614 0.0892532 0.0887277 0.0901797 0.0905337 0.0899614 0.0906789 0.0917796 0.0930459 0.0947894 0.0976045 0.0985857 0.1003443 0.1007281 0.1009501 0.1018683 0.1032068 0.1047455 0.1065945 0.1075091 0.1078353 0.1061196 0.1065857 0.1087148 0.1089483 0.1098117 0.1105224 0.1095223 0.1082381 0.1067355 0.1069765 0.1071878 0.1082369 0.1093421 0.109819 0.11185 0.1144274 0.1167922 0.1180822 0.1186669 0.1204462 0.1198099 0.1185896 0.1132963 0.1118587 0.1091871 0.1050476 0.102304 0.10089 0.1004672 0.1011395 0.1060796 0.1129575 0.1175725 0.1232193 0.130118 0.135125 0.1360825 0.1359793 0.1326176 0.1321644 0.129457 0.124572 0.1217336 0.1186226 0.1168418 0.1160937 0.1179839 0.1201002 0.1234569 0.1275045 0.1299747 0.1325541 0.1335283 0.1330306 0.1334459 0.1314906 0.1329787 0.1345618 0.1338015 0.1310492 0.1291785 0.1258195 0.122067 0.1202742 0.1186862 0.1175583 0.1190172 0.1204426 0.1243681 0.1280569 0.1351692 0.140639 0.1445625 0.1505314 0.1560489 0.1616585 0.1692669 0.1741895 0.1770809 0.1789952 0.1818543 0.1759706 0.1759706 0.16862 0.1612694 0.1539188 0.1465682 0.1392176 0.131867 0.1245163 0.1171657 0.1098151 0.1024645 0.0951139 0.0877633 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0805511 0.0806895 0.0808278 0.080966 0.081104 0.0812418 0.0813793 0.0815166 0.0816535 0.0817901 0.0819263 0.082062 0.0821972 0.0823319 0.0824661 0.0825996 0.0827325 0.0828646 0.0829961 0.0831268 0.0832567 0.0833857 0.0835139 0.0836412 0.0837675 0.0838928 0.084017 0.0841402 0.0842623 0.0843833 0.0845031 0.0846216 0.0847389 0.084855 0.0849697 0.0850831 0.0851951 0.0853056 0.0854147 0.0855224 0.0856285 0.0857331 0.0858361 0.0859375 0.0860372 0.0861353 0.0862317 0.0863264 0.0864193 0.0865104 0.0865998 0.0866873 0.0867729 0.0868567 0.0869385 0.0870184 0.0870964 0.0871723 0.0872463 0.0873182 0.0873881 0.087456 0.0875217 0.0875854 0.0876469 0.0877062 0.0877634 0.0878185 0.0878713 0.0879219 0.0879703 0.0880165 0.0880604 0.088102 0.0881414 0.0881785 0.0882132 0.0882457 0.0882758 0.0883036 0.0883291 0.0883522 0.088373 0.0883914 0.0884075 0.0884211 0.0884325 0.0884414 0.0884479 0.0884521 0.0884539 0.0884533 0.0884503 0.088445 0.0884372 0.0884271 0.0884146 0.0883997 0.0883825 0.0883629 0.088341 0.0883167 0.08829 0.0882611 0.0882298 0.0881961 0.0881602 0.088122 0.0880815 0.0880387 0.0879937 0.0879464 0.0878969 0.0878452 0.0877912 0.0877351 0.0876768 0.0876164 0.0875538 0.0874891 0.0874223 0.0873535 0.0872825 0.0872096 0.0871346 0.0870576 0.0869787 0.0868978 0.086815 0.0867303 0.0866437 0.0865553 0.0864651 0.0863731 0.0862793 0.0861837 0.0860865 0.0859876 0.085887 0.0857848 0.085681 0.0855756 0.0854687 0.0853604 0.0852505 0.0851392 0.0850266 0.0849125 0.0847971 0.0846804 0.0845625 0.0844433 0.0843229 0.0842014 0.0840788 0.083955 0.0838302 0.0837044 0.0835777 0.0834499 0.0833213 0.0831919 0.0830616 0.0829305 0.0827986 0.0826661 0.0825329 0.0823991 0.0822646 0.0821297 0.0819942 0.0818582 0.0817219 0.0815851 0.081448 0.0813106 0.0811729 0.081035 0.0808969 0.0807586 0.0806203 0.0804819 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0804127 0.0805511 0.0806895 0.0808278 0.080966 0.081104 0.0812418 0.0813793 0.0815166 0.0816535 0.0817901 0.0819263 0.082062 0.0821972 0.0823319 0.0824661 0.0825996 0.0827325 0.0828646 0.0829961 0.0831268 0.0832567 0.0833857 0.0835139 0.0836412 0.0837675 0.0838928 0.084017 0.0841402 0.0842623 0.0843833 0.0845031 0.0846216 0.0847389 0.084855 0.0849697 0.0850831 0.0851951 0.0853056 0.0854147 0.0855224 0.0856285 0.0857331 0.0858361 0.0859375 0.0860372 0.0861353 0.0862317 0.0863264 0.0864193 0.0865104 0.0865998 0.0866873 0.0867729 0.0868567 0.0869385 0.0870184 0.0870964 0.0871723 0.0872463 0.0873182 0.0873881 0.087456 0.0875217 0.0875854 0.0876469 0.0877062 0.0877634 0.0878185 0.0878713 0.0879219 0.0879703 0.0880165 0.0880604 0.088102 0.0881414 0.0881785 0.0882132"


## vaccination stratgies:
##      0. no vaccine
##      1. random: "-vac1-frac-10 0.11111111 -vac1-frac-20 0.11111111 -vac1-frac-30 0.11111111 -vac1-frac-40 0.11111111 -vac1-frac-50 0.11111111 -vac1-frac-60 0.11111111 -vac1-frac-70 0.11111111 -vac1-frac-80 0.11111111" 
##      2. 16-29 only*: "-vac1-frac-10 0.28571400 -vac1-frac-20 0.71428500 -vac1-frac-30 0.00000014 -vac1-frac-40 0.00000014 -vac1-frac-50 0.00000014 -vac1-frac-60 0.00000014 -vac1-frac-70 0.00000014 -vac1-frac-80 0.00000014"
##      3. 30-59 only*: "-vac1-frac-10 0.00000017 -vac1-frac-20 0.00000017 -vac1-frac-30 0.33333300 -vac1-frac-40 0.33333300 -vac1-frac-50 0.33333300 -vac1-frac-60 0.00000017 -vac1-frac-70 0.00000017 -vac1-frac-80 0.00000017"
##      4. 60+ only*: "-vac1-frac-10 0.00000017 -vac1-frac-20 0.00000017 -vac1-frac-30 0.00000017 -vac1-frac-40 0.00000017 -vac1-frac-50 0.00000017 -vac1-frac-60 0.33333300 -vac1-frac-70 0.33333300 -vac1-frac-80 0.33333300"
##      5. random 50/50 excl. 30-59*: "-vac1-frac-10 0.16666650 -vac1-frac-20 0.16666650 -vac1-frac-30 0.00000033 -vac1-frac-40 0.00000033 -vac1-frac-50 0.00000033 -vac1-frac-60 0.16666650 -vac1-frac-70 0.16666650 -vac1-frac-80 0.16666650"
##      6. random 25/75 excl. 30-59*: "-vac1-frac-10 0.08333325 -vac1-frac-20 0.08333325 -vac1-frac-30 0.00000033 -vac1-frac-40 0.00000033 -vac1-frac-50 0.00000033 -vac1-frac-60 0.24999975 -vac1-frac-70 0.24999975 -vac1-frac-80 0.24999975"
##      7. random 75/25 excl. 30-59*: "-vac1-frac-10 0.24999975 -vac1-frac-20 0.24999975 -vac1-frac-30 0.00000033 -vac1-frac-40 0.00000033 -vac1-frac-50 0.00000033 -vac1-frac-60 0.08333325 -vac1-frac-70 0.08333325 -vac1-frac-80 0.08333325"
##      8. random among 20-49 and 70+*: "-vac1-frac-10 0.00000025 -vac1-frac-20 0.1999998 -vac1-frac-30 0.1999998 -vac1-frac-40 0.1999998 -vac1-frac-50 0.00000025 -vac1-frac-60 0.00000025 -vac1-frac-70 0.1999998 -vac1-frac-80 0.1999998"
##      9. random 50/50 among 20-49 and 70+*: "-vac1-frac-10 0.00000025 -vac1-frac-20 0.16666650 -vac1-frac-30 0.16666650 -vac1-frac-40 0.16666650 -vac1-frac-50 0.00000025 -vac1-frac-60 0.00000025 -vac1-frac-70 0.24999975 -vac1-frac-80 0.24999975"
##     10. random 25/75 among 20-49 and 70+*: "-vac1-frac-10 0.00000025 -vac1-frac-20 0.08333325 -vac1-frac-30 0.08333325 -vac1-frac-40 0.08333325 -vac1-frac-50 0.00000025 -vac1-frac-60 0.00000025 -vac1-frac-70 0.374999625 -vac1-frac-80 0.374999625"
##     11. random 75/25 among 20-49 and 70+*: "-vac1-frac-10 0.00000025 -vac1-frac-20 0.24999975 -vac1-frac-30 0.24999975 -vac1-frac-40 0.24999975 -vac1-frac-50 0.00000025 -vac1-frac-60 0.00000025 -vac1-frac-70 0.124999875 -vac1-frac-80 0.124999875"
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