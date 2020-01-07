#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THSD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THSD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THSD7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THSD7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THTPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THUMPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THUMPD2 
