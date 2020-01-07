#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=30gb
#SBATCH -p standard
#SBATCH -A zanglab

#SBATCH --mail-user=zhenjia@virginia.edu
#SBATCH --mail-type=end
#SBATCH --job-name=HiCpro_s1_prostate_Bglii
#SBATCH --export=ALL
#SBATCH --array=1-11

FASTQFILE=$SLURM_SUBMIT_DIR/inputfiles_prostate_Bglii.txt; export FASTQFILE
make --file /nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f0_HiC_process/prostate_HiC_Pro_BglII/config_prostate_Bglii.txt CONFIG_SYS=/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/config-system.txt all_sub 2>&1
