#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NELL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NELL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEMF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NENF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEO1 
