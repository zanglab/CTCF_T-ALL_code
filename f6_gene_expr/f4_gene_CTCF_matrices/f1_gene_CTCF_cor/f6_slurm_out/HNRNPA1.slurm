#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HNRNPA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA1P10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA1P33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA2B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA3 
