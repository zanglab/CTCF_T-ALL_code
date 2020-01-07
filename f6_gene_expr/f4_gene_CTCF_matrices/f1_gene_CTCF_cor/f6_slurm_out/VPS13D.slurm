#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VPS13D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS13D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS26A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS26B 
