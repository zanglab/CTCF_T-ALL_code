#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LEF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEFTY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEFTY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEKR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LELP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEMD1 
