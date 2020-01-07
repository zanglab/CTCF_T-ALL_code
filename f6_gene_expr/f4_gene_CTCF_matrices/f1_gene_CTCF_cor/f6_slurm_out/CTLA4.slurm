#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTLA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTLA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTNNB1 
