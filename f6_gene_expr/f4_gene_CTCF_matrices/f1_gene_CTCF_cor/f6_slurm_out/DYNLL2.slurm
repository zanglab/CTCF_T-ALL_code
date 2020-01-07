#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DYNLL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNLL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNLRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNLRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNLT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNLT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYRK1A 
