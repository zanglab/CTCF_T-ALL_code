#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRA2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRA2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRABD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRABD2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRABD2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRADD 
