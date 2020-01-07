#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LPIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LPIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPP 
