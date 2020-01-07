#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSPB7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPBAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPD1 
