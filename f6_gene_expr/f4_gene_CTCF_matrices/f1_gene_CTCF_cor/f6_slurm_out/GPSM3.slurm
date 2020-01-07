#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPSM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPSM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX3 
