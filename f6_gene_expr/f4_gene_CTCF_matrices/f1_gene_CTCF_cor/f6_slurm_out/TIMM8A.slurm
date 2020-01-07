#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIMM8A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMMDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMP2 
