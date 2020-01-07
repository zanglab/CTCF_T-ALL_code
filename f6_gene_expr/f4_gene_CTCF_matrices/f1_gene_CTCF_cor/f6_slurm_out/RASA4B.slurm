#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RASA4B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RASA4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASA4CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASD1 
