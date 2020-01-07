#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C9orf24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf50 
