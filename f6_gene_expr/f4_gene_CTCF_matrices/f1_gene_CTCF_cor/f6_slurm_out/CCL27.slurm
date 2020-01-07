#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCL27.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL3L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL4 
