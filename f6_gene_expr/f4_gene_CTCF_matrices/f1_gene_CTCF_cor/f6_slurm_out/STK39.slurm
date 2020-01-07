#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STK39.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STK39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STKLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STMN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STMN2 
