#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPALPP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPALPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPANK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPAT4 
