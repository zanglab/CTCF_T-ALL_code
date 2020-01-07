#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAI14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAI2 
