#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TFAP2D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAP2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAP2E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFB1M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFB2M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFCP2 
