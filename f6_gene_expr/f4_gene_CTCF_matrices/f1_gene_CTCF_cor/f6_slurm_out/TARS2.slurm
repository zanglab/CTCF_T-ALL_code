#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TARS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TARSL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS1R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS1R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS1R3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TAS2R1 
