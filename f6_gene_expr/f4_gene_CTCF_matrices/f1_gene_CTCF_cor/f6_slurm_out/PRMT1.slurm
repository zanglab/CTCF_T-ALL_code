#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRMT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRMT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRMT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRMT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRMT7 
