#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NARF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NARF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NARFL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NASP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAT1 
