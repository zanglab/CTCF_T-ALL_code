#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HNRNPK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPKP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPLL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPR 
