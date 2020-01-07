#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATA25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31A5 
