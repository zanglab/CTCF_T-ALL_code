#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEK9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NELFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NELFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NELFCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NELFE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NELL1 
