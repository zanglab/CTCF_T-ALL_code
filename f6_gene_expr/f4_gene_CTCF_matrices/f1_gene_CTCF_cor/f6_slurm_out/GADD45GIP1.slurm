#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GADD45GIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GADD45GIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GADL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE12B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE12D 
