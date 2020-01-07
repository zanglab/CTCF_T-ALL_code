#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAS2R60.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAT 
