#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAGA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAGK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAGLU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAGPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAGS 
