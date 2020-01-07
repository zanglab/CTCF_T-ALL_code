#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCBP3 
