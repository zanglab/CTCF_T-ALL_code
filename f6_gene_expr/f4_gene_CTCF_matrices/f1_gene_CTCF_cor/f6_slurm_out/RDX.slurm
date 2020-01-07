#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RDX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RDX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REC114 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RECK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RECQL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RECQL4 
