#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THPO.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THRAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THRSP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THSD1 
