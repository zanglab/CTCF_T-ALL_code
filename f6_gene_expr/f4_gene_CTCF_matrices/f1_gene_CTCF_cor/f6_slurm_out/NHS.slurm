#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NHS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NHS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHSL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NICN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NID1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NID2 
