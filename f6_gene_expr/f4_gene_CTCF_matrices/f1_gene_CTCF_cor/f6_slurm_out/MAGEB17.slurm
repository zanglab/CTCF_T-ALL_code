#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGEB17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEB6 
