#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RARS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASA4 
