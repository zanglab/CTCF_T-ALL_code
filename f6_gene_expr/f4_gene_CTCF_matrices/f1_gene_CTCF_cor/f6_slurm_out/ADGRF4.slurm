#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADGRF4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG4 
