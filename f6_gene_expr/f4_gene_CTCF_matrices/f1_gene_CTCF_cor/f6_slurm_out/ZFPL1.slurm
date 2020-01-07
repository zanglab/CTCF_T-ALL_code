#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFPL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFPM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFPM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFX 
