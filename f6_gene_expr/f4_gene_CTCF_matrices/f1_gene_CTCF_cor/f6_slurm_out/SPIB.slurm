#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPIB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPICE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIDR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIN2A 
