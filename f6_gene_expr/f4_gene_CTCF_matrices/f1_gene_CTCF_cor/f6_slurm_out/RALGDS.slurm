#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RALGDS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RALGDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALGPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALGPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALYL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAMP1 
