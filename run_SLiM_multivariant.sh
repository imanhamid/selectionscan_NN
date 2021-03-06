#!/bin/bash
#SBATCH -p scavenger
#SBATCH --job-name="multivariant"
#SBATCH -a 1-100

c=$SLURM_ARRAY_TASK_ID

declare -i c10="${c}0"

for i in {0..9};
do

c_array[$i]=$((c10 - i))

done

for i in "${c_array[@]}"
do

~/home/SLiM_build/slim -d out=\'"/work/ih49/simulations/segmentation/multivariant/multivariant_seed-$i"\' ~/home/selectionscan_NN/multivariant.slim

done
