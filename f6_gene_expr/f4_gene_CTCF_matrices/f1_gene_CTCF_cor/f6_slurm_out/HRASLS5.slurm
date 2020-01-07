#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HRASLS5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HRASLS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HRC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HRCT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HRG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HRH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HRH2 
