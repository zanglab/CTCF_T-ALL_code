#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DIMT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DIMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIP2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIP2B 
