#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KANK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KANK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANSL1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANSL2 
