#!/bin/bash


loc="RI"
declare -a dir_pre_a=("20201226-med-allvac" "20201226-low-allvac" "20201226-high-allvac")
# declare -a dir_pre_a=("20201226-med-sens204970" "20201226-low-sens204970" "20201226-high-sens204970")

declare -a dir_sfx=("tot_50k" "tot_50k-notest" "tot_300k" "tot_300k-notest")
# declare -a dir_sfx=("tot_50k" "tot_300k")


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
