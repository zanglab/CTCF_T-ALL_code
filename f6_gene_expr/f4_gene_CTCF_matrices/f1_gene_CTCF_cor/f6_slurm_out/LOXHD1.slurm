#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LOXHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LOXHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LOXL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LOXL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LOXL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LOXL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPA 
