#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1orf68.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf94 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf141 
