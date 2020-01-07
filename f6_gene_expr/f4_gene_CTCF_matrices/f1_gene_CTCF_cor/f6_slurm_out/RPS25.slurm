#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS26P11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS27A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS27L 
