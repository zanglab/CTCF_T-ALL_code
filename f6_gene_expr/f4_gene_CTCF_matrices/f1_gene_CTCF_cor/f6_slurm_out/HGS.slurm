#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HGS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HGS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HGSNAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHATL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHEX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHIP 
