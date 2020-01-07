#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UMAD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UMAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UMOD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UMODL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UMPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC119 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC119B 
