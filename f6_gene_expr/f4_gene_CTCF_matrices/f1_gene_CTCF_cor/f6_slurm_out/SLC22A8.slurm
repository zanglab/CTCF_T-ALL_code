#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC22A8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC23A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC23A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC23A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC24A1 
