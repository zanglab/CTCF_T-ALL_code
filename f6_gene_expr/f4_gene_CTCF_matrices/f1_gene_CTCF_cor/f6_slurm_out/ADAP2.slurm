#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADARB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADARB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAT2 
