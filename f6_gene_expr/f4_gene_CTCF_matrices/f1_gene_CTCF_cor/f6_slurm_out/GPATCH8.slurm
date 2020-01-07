#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPATCH8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPBAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPBP1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPC2 
