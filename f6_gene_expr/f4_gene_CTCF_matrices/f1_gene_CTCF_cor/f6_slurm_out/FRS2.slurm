#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FRS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FRS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRYL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRZB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSBP 
