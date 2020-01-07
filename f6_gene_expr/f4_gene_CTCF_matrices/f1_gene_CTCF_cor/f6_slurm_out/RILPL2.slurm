#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RILPL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RILPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMBP3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMBP3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIMKLA 
