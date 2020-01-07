#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHGDH 
