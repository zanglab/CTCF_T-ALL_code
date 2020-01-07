#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RREB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RREB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRM2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRN3 
