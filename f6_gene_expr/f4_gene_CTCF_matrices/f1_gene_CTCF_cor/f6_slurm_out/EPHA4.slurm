#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPHA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHB1 
