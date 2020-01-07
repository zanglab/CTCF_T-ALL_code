#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRKAR2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAR2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAR2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCE 
