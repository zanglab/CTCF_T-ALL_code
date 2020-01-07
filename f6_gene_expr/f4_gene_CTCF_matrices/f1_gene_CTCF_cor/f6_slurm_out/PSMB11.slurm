#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMB11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMB6 
