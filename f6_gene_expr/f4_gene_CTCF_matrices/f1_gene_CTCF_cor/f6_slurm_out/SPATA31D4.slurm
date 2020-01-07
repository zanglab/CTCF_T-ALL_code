#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATA31D4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31D4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31D5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA4 
