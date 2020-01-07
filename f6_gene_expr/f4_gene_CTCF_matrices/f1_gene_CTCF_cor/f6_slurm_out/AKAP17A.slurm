#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKAP17A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP17A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKAP6 
