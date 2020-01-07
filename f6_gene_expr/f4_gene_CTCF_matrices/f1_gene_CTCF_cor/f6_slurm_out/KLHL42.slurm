#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHL42.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL9 
