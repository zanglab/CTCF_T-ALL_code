#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIGM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGR 
