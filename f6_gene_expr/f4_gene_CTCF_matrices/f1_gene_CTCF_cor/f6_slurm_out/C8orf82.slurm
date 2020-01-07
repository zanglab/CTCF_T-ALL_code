#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C8orf82.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf86 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf89 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C9orf116 
