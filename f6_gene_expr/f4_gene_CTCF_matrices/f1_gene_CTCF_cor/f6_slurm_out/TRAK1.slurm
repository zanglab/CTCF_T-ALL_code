#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRAK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAM1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRANK1 
