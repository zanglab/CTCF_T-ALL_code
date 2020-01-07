#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SQLE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SQSTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRBD1 
