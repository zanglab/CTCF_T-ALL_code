#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFHX3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFHX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFHX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFP28 
