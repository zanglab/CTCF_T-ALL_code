#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SAT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SATB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SATB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SATL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAXO1 
