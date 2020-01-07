#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RP9P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPA3 
