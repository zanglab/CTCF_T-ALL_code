#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACOX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOXL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACP2 
