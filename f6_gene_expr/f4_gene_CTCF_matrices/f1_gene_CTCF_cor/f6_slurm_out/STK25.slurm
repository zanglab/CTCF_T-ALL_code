#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STK25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STK25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK32A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK32B 
