#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NANOS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NANOS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NANOS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NANOS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NANP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NANS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAP1L1 
