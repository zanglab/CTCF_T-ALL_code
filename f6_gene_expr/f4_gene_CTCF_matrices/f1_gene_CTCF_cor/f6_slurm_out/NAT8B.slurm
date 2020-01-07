#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAT8B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT8L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NATD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAV2 
