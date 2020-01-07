#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NSMCE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NSMCE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSMCE4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSMF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN3 
