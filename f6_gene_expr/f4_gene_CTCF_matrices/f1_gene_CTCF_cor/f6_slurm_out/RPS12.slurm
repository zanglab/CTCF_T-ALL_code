#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS15A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS15AP10 
