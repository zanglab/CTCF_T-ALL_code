#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFS6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFV2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFV3 
