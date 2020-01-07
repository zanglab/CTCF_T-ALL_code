#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FCGRT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCHO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCHO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCHSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCHSD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCMR 
