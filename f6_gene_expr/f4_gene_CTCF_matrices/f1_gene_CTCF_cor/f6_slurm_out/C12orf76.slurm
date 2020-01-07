#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C12orf76.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf119 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf132 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf177 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf178 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C14orf180 
