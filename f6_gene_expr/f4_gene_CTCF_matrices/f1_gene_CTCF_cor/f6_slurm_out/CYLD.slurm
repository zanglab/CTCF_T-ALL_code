#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYLD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYLD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP11A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP11B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP11B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP17A1 
