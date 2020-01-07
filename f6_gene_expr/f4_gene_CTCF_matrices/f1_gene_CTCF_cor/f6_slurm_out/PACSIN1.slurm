#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PACSIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PACSIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PACSIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PACSIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PADI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PADI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PADI3 
