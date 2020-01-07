#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAS2R10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R3 
