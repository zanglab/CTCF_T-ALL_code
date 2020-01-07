#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNG8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNGT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNL3 
