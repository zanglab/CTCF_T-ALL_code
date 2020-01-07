#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5AR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5AR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf24 
