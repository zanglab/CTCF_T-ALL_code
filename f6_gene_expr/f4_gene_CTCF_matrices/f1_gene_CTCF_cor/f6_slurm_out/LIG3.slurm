#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LIG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LIG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRA4 
