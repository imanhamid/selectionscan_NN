#!/bin/bash
#SBATCH -p scavenger
#SBATCH --job-name="png"
#SBATCH -a 1-10
#SBATCH --mem=5G
#SBATCH -o "100_png.out"

module load R/4.0.0

c=$SLURM_ARRAY_TASK_ID


declare -i c100="${c}00"

for i in {0..99};
do

c_array[$i]=$((c100 - i))

done

for i in "${c_array[@]}"
do


filename=$(ls /work/ih49/simulations/segmentation/two-strong/*_seed-${i}_alltracts.txt | head -n 1)
~/home/selectionscan_NN/2strong_segmentation.R $filename

done
