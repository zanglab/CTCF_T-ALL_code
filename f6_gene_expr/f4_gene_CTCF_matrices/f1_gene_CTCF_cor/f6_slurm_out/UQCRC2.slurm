#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UQCRC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCRC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCRFS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCRH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCRHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UQCRQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i URAD 
