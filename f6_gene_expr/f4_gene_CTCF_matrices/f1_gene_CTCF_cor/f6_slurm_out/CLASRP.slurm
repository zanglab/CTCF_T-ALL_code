#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLASRP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLASRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCA3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCA4 
