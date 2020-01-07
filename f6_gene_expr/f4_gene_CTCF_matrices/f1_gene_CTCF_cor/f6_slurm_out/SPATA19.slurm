#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATA19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA24 
