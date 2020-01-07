#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NTSR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NTSR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUAK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUAK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUBP2 
