#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C14orf28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf93 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf38-AP3S2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf39 
