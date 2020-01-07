#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR10Z1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10Z1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11H4 
