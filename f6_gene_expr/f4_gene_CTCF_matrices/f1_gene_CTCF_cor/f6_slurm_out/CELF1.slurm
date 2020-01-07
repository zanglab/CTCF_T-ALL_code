#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CELF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CELF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELF6 
