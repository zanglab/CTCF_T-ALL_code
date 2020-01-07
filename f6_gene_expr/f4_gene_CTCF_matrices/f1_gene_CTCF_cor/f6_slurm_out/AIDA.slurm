#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AIDA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AIDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIF1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIFM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIFM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AIFM3 
