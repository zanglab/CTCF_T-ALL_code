#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTTN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTTNBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTTNBP2NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTU2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTXN1 
