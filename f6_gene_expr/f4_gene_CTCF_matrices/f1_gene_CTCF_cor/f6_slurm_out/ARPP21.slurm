#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARPP21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPP21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRDC2 
