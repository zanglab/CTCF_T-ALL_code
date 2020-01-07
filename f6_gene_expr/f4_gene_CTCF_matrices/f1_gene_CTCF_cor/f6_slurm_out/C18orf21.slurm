#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C18orf21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C18orf21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C18orf25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C18orf32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C18orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C18orf63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C19orf12 
