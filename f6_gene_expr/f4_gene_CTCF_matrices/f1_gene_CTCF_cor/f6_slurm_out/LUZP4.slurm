#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LUZP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LUZP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LUZP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LVRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LXN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6E 
