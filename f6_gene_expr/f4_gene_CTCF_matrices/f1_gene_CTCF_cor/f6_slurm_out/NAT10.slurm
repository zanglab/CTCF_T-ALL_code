#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAT10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT8 
