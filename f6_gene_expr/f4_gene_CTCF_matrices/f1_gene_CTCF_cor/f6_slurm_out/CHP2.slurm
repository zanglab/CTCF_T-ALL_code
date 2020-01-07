#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHPF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRD 
