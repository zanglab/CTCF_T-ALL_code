#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCDC2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCDC2C 
