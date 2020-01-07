#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C16orf87.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf89 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf90 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf92 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf95 
