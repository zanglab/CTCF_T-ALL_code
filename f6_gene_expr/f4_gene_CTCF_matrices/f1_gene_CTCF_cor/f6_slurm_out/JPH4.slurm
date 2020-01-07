#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o JPH4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i JPH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JRK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JRKL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JSRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JUN 
