#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGEF10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF10L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF16 
