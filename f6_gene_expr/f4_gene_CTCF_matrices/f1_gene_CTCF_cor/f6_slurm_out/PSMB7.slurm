#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMB7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMC3 
