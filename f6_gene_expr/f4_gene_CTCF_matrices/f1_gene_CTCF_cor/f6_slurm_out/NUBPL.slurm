#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUBPL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUBPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUCB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUCB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUCKS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDCD1 
