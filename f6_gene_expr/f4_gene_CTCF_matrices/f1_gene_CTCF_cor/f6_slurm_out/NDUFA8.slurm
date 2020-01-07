#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFA8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFAF3 
