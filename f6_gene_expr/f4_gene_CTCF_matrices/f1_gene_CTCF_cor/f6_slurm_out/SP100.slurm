#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SP100.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SP100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP110 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP140 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP140L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP3 
