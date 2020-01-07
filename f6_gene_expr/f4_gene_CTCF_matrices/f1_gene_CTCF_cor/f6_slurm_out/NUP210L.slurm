#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUP210L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP210L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP210P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP214 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP43 
