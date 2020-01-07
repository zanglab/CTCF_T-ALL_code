#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMNDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMNDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMOC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMOC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPD1 
