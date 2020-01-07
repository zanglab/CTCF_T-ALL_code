#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRCP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRCP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRCT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB3L2 
