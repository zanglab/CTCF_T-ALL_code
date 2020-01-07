#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDZD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD9 
