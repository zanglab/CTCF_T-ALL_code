#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD18A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD18B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD18DP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD19P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD2 
