#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o URAHP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i URAHP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i URB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i URB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i URGCP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i URGCP-MRPS24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i URI1 
