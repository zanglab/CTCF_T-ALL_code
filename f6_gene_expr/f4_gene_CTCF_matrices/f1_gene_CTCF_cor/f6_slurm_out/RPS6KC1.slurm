#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS6KC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS7P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS9 
