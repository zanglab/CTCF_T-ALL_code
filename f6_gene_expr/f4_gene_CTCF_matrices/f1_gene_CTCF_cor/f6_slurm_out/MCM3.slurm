#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM3AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM7 
