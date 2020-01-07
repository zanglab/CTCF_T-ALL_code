#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OPRK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OPRK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPRL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPTC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10A2 
