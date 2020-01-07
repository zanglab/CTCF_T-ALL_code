#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASAH2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASAH2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASB10 
