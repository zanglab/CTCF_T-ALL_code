#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GKAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GKAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GKN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GKN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLB1L 
