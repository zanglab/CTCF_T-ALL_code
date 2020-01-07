#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGAP25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP30 
