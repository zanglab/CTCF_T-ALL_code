#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ULK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ULK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ULK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ULK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ULK4P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ULK4P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ULK4P3 
