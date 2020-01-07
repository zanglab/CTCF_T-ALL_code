#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFITM5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFITM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNA22P 
