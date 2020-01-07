#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TGS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TGS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THADA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THAP11 
