#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APOPT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APOPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APPBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APRT 
