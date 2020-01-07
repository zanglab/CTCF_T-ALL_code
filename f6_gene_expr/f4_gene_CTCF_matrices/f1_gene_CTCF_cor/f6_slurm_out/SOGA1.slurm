#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOGA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOGA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOHLH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOHLH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SON 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORBS1 
