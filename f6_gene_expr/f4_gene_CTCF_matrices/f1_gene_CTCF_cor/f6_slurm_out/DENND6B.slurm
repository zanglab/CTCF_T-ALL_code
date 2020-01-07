#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DENND6B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEPDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEPDC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEPDC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEPDC5 
