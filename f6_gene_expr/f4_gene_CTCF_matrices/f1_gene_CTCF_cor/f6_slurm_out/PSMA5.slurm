#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMA5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB10 
