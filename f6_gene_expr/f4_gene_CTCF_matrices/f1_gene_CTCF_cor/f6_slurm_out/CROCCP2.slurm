#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CROCCP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CROCCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CROCCP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CROT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRTAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRTAM 
