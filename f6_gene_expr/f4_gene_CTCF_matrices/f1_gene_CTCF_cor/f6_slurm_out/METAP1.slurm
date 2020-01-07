#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o METAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i METAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METAP1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METRNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL1 
