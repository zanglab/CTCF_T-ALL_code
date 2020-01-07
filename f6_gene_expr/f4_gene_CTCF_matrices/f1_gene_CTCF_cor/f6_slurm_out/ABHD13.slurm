#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABHD13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD14A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD14A-ACY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD14B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD16A 
