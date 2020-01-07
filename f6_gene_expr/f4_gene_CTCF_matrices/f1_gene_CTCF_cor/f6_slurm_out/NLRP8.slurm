#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NLRP8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMBR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMD3 
