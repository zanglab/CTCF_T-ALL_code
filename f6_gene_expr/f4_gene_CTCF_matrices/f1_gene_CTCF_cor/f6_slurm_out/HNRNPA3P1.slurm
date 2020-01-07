#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HNRNPA3P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPA3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPCL4 
