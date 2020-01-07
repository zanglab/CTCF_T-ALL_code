#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o YBX3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i YBX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YBX3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YDJC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YEATS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YEATS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YES1 
