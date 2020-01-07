#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRCH4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRCH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRCOL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRFN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRFN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRFN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRFN4 
