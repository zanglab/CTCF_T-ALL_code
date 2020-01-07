#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CREB3L3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB3L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB3L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREBBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREBL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CREBRF 
