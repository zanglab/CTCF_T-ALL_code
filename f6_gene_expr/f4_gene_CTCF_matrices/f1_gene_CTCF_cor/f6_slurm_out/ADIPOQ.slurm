#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADIPOQ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADIPOQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADIPOR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADIPOR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADIRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADM 
