#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SGIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SGIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGK494 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGMS1 
