#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DLK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DLK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLSTP1 
