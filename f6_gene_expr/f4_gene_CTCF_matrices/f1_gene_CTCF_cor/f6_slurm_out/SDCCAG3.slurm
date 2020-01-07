#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SDCCAG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SDCCAG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDCCAG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDF2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDF4 
