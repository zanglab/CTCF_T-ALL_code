#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATA31A6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA31D3 
