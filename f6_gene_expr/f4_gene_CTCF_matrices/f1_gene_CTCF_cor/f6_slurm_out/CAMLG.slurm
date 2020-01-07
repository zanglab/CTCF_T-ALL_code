#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAMLG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMSAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMSAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMSAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAMTA1 
