#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SHB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHBG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHC4 
