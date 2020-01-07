#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCHE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCHE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCKDHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCKDHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCKDK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL11A 
