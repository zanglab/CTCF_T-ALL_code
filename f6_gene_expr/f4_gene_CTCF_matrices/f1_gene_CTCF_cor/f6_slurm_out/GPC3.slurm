#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPCPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPD1 
