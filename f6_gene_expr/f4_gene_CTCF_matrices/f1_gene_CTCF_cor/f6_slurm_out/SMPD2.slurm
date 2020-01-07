#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMPD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPDL3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPDL3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMPX 
