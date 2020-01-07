#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NRG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NRG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRIP3 
