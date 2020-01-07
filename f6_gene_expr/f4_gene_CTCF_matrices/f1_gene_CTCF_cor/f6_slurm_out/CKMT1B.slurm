#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CKMT1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CKMT1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKS1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLASP2 
