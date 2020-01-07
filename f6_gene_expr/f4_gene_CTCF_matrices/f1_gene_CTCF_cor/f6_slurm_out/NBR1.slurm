#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NBR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NBR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCALD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAPD2 
