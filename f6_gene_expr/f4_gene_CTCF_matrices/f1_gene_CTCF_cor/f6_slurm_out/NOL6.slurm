#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOL6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOM1 
