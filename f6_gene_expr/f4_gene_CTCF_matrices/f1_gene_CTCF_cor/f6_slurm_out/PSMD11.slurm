#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMD11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD3 
