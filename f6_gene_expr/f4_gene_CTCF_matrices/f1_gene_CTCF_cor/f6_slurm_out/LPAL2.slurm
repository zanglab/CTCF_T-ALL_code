#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LPAL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAR5 
