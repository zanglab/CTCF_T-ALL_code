#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TM4SF19-TCTEX1D2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF19-TCTEX1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM6SF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM6SF2 
