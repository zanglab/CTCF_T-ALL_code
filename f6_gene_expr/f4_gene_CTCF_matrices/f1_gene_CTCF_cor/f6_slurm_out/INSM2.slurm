#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INSM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INSM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSRR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS12 
