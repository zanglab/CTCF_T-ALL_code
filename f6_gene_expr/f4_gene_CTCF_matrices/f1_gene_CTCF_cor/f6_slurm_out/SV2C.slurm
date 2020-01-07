#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SV2C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SV2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVIL-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SVILP1 
