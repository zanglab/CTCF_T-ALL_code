#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASH2L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASH2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASIC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASIC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASIC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASIC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASIC5 
