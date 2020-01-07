#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARNTL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARNTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARNTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC3 
