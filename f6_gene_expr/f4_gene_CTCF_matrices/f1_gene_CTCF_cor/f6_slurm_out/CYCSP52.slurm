#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYCSP52.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYCSP52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYFIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYFIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYHR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYLC2 
