#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPM1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1H 
