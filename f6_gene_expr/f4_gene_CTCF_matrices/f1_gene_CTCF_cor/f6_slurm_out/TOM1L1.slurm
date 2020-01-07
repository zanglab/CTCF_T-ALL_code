#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOM1L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOM1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOM1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM20L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOMM34 
