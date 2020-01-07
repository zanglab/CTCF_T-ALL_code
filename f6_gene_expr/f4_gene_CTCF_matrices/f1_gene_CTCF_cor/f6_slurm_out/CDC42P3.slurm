#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDC42P3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42SE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42SE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC6 
