#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOL10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL4L 
