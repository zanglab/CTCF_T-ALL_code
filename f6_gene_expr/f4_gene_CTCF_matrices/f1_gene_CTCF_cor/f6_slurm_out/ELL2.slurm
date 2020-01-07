#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMOD1 
