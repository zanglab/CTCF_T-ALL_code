#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FER1L4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FER1L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FER1L5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FER1L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FERD3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FERMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FERMT2 
