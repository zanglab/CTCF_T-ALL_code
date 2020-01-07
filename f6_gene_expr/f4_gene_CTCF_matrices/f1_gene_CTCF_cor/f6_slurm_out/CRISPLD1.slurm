#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRISPLD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRISPLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRISPLD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRKL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRLF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRLF2 
