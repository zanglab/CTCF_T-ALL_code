#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEC24C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC24C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC24D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC31A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC31B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC61A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC61A2 
