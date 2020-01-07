#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACY3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACYP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACYP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAD2 
