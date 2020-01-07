#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZEB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZEB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND3 
