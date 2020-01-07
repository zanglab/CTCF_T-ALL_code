#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ROMO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ROMO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROPN1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROPN1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROR2 
