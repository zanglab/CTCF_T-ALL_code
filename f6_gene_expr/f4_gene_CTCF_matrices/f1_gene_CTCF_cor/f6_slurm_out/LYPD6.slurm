#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LYPD6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPLA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPLA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPLA2P2 
