#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAPB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPEPLD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPSA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPSB 
