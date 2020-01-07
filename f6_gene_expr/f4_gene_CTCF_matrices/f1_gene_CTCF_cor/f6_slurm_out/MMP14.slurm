#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MMP14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP2 
