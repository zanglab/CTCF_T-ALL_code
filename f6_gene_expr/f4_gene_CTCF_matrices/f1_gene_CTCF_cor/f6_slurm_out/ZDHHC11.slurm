#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZDHHC11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC16 
