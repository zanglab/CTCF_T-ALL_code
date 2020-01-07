#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE26 
