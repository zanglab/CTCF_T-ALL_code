#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NRK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NRK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRN1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRP1 
