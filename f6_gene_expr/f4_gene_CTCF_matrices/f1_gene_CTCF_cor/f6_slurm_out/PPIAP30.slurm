#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPIAP30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIAP30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPID 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIEL 
