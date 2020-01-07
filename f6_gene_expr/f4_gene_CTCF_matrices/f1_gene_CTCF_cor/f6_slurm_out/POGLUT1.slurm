#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POGLUT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POGLUT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POGZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLD1 
