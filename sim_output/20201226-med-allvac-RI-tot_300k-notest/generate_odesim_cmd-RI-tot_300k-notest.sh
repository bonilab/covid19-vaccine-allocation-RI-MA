#!/bin/bash

##
## script to generate a series of run commands
## and create empty output folders to hold run outputs from those commands
## all commands will be save to cmd.sh
##

script="./run_vac_strategies-RI-notest.sh"
odesim_dir="../../cpp-v6-test-vaccination/"

declare -a nat_imm_a=(540)
declare -a normalcyday_a=(1098)
declare -a vac_duration_a=(540)
declare -a vac_halflife_a=(180 180 180 360 360 540 540)
declare -a vac_slope_a=(2 3 4 1.5 2 1.5 2)


startday="370 400"
endday="399 429"
dpd="5000 5000"

sd="\""$startday"\""
ed="\""$endday"\""
dp="\""$dpd"\""

exec 3<> cmd.sh
echo "#!/bin/bash" >&3
for nat_imm in "${nat_imm_a[@]}"
do
    for normalcyday in "${normalcyday_a[@]}"
    do
        for vac_duration in "${vac_duration_a[@]}"
        do
            i=0
            for vac_halflife in "${vac_halflife_a[@]}"
            do
                vac_slope=${vac_slope_a[$i]}

                outdir="s${startday:0:3}_e${endday:((${#endday}-3)):3}_d${dpd:0:4}_i${nat_imm}_c${normalcyday}_u${vac_duration}_t${vac_halflife}_p${vac_slope}"
                if [ -d "$outdir" ]; then
                    rm $outdir/*
                else
                    mkdir $outdir
                fi
                echo "" >&3
                echo "$script" "-r" $odesim_dir "-o" $outdir "-s" "$sd" "-e" "$ed" "-d" "$dp" "-i" $nat_imm "-c" $normalcyday "-u" $vac_duration "-t" $vac_halflife "-p" $vac_slope >&3

                i=$(($i+1))
            done
        done
    done
done
exec 3>&-

exit
