#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALPL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALPPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALS2CL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALS2CR12 
