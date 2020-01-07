#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RASL10A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RASL10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASL10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASL11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASL11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASL12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF1 
