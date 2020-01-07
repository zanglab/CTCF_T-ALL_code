#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARV1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARVCF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AS3MT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASAH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASAH2 
