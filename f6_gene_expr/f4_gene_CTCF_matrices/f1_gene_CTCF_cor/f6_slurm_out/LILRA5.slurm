#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LILRA5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRB4 
