#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOMM40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM40L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TONSL 
