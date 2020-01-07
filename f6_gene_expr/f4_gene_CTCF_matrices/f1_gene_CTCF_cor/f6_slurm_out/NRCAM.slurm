#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NRCAM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NRCAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRDE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NREP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRG2 
