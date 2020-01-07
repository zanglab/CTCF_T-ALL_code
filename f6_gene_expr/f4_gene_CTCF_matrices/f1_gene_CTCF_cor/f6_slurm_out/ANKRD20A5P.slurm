#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD20A5P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD20A5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD20A8P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD20A9P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD24 
