#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAGLN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAGLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAGLN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TALDO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAMM41 
