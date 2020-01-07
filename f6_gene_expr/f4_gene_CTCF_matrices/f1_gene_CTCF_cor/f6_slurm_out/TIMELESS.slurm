#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIMELESS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMELESS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM17A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMM17B 
