#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMR3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMR3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMTNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMTNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMU1 
