#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ECHDC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ECHDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECHS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECM2 
