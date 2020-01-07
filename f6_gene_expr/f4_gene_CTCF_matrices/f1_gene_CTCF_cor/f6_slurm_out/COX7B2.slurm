#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COX7B2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COX7B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX7C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX8C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPA1 
