#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DHX33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX38 
