#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FANCA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCD2OS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCE 
