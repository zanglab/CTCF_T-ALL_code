#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRSF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF7 
