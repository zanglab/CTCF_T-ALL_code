#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MOK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MOK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MON1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MON1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MON2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORC2 
