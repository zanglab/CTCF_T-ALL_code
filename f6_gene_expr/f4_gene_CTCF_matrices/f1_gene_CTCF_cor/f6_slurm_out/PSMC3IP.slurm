#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMC3IP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC3IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD10 
