#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HPDL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HPDL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPGD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPGDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPRT1 
