#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AP1G2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1M2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1S2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1S3 
