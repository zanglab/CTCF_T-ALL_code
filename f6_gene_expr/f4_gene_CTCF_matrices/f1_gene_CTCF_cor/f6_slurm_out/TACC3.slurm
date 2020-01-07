#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TACC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TACC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TACSTD2 
