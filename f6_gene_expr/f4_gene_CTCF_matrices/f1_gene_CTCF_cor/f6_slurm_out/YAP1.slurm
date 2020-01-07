#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o YAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i YAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YBEY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YBX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YBX2 
