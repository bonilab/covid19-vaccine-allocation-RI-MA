#!/bin/bash

loc="MA"

## similar to Figure 3, 5, 7 but for MA
declare -a dir_pre_a=("20210109-med-allvac" "20210109-low-allvac" "20210109-high-allvac")
declare -a dir_sfx=("tot_300k" "tot_300k-notest" "tot_1800k" "tot_1800k-notest")

## similar to Figure 6 but for MA
# declare -a dir_pre_a=("20210109-med-sens204970" "20210109-low-sens204970" "20210109-high-sens204970")
# declare -a dir_sfx=("tot_300k" "tot_1800k")

for dir_pre in "${dir_pre_a[@]}"
do
	for ds in "${dir_sfx[@]}"
	do
		lds="${loc}-${ds}"
		dir_name="${dir_pre}-${lds}"
		gen_scrpt="generate_odesim_cmd-${lds}.sh"
		
		echo ""
		echo "directory $dir_name"
		cd $dir_name
		
		rm *.dat

		if [ cmd.sh ]; then
		rm cmd.sh
		fi
		echo "generate cmd"
		./$gen_scrpt
		chmod +x cmd.sh
		

		echo "run cmd"
		./cmd.sh

		cd ..
	done
done

exit
