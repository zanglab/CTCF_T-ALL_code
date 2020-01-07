#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NRP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRROS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRSN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRSN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRXN1 
