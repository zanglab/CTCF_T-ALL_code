#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IREB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IREB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF2BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF2BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF2BPL 
