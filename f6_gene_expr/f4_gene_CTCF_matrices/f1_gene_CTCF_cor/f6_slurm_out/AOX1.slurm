#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AOX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AOX2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1AR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1B1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP1G1 
