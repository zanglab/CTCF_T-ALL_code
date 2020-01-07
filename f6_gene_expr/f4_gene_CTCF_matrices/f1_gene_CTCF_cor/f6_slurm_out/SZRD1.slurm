#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SZRD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SZRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SZT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAAR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAB2 
