#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAP2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFAP2C 
