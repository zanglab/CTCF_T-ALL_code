#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THOC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THOC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THOC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THOC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THOC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THOC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THOP1 
