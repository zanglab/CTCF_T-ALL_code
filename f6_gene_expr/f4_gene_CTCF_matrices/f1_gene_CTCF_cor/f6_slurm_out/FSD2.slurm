#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FSD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FSD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FST 
