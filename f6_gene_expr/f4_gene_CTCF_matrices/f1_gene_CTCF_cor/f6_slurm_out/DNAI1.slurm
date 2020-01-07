#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJA1P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJA3 
