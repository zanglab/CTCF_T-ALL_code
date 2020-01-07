#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ST3GAL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ST3GAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST3GAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST3GAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST3GAL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST3GAL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST3GAL6 
