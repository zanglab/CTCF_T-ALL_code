#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EEFSEC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EEFSEC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB12 
