#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o py1_slurms/slurm_chr3.out

time python 1_get_completed_view2M_interaction.py -c chr3 
