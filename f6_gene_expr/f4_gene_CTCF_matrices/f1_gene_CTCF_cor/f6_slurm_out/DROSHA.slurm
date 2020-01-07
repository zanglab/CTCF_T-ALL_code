#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DROSHA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DROSHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSCAM 
