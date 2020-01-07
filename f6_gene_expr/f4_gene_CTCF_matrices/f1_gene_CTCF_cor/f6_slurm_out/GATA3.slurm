#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GATA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GATA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATAD2A 
