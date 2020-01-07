#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS16P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS18P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS19 
