#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTNNBIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNBIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNND2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTPS1 
