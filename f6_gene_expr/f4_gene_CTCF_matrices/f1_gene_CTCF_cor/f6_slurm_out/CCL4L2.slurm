#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCL4L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL4L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCM2L 
