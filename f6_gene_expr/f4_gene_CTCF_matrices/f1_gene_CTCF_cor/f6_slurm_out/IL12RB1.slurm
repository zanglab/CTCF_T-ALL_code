#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL12RB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL12RB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL12RB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL13RA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL13RA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL15 
