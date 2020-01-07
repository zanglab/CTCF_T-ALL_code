#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZSCAN12P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN12P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN21 
