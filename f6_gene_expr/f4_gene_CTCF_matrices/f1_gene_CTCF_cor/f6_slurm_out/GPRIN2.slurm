#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPRIN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPSM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPSM2 
