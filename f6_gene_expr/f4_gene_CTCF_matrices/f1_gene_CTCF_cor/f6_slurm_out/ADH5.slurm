#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADH5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADHFE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADIG 
