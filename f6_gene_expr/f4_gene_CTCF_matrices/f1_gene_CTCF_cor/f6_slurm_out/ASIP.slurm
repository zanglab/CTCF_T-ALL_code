#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASMTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASNS 
