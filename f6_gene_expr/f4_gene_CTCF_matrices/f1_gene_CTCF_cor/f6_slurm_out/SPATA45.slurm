#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATA45.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA5L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA6L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATA7 
