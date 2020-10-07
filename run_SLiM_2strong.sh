#!/bin/bash
#SBATCH -p scavenger
#SBATCH --job-name="2strong"
#SBATCH -a 1-10

c=$SLURM_ARRAY_TASK_ID

declare -i c100="${c}00"

for i in {0..99};
do

c_array[$i]=$((c100 - i))

done

for i in "${c_array[@]}"
do

~/home/SLiM_build/slim -d out=\'"/work/ih49/simulations/segmentation/two-strong/two-strong_seed-$i"\' ~/home/selectionscan_NN/2strong.slim

done
