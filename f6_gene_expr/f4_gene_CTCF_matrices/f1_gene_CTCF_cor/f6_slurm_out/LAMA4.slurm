#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LAMA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMB4 
