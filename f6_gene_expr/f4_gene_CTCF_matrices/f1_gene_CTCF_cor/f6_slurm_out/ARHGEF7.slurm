#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGEF7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID3A 
