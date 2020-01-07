#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MSL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MSL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSL3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSMO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSMP 
