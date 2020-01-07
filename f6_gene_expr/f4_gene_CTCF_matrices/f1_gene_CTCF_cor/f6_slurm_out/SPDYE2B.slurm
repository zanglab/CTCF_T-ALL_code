#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPDYE2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE7P 
