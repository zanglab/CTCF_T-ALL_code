#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RABEPK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RABEPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABGAP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABGEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABGGTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABGGTB 
