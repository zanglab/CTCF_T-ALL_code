#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDZK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZK1IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZRN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZRN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEA15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEAK1 
