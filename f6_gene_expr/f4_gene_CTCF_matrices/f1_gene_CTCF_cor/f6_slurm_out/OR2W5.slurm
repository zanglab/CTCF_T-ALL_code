#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2W5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2W5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2Z1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR3A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR3A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR3A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR3A4P 
