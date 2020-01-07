#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RTN4IP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN4IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN4R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN4RL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN4RL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTP2 
