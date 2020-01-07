#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMY1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMY1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMY1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMY2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMY2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMZ2 
