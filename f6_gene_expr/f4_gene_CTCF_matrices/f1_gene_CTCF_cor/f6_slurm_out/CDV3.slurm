#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDV3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDV3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDY2A 
