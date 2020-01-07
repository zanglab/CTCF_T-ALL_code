#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSPB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB2-C11orf52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB6 
