#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LUC7L3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LUC7L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LUM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LURAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LURAP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LUZP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LUZP2 
