#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLB1L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLB1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLB1L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLCCI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLDN 
