#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MLX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MLX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLXIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLXIPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLYCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMAA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMAB 
