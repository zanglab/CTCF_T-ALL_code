#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPSA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPSA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPSAP52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPSAP58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPSAP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPTOR 
