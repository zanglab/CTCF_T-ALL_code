#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OTUB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTUD5 
