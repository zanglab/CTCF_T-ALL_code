#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NRAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NRAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRARP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRBF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRBP2 
