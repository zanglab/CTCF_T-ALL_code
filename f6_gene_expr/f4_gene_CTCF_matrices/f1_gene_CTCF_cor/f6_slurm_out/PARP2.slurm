#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PARP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PARP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARP9 
