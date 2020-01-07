#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STK32C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STK32C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK38L 
