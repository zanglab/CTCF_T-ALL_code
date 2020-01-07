#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NTPCR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NTPCR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTRK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTRK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTRK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTSR1 
