#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCTS2P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCTS2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCUR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDFI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MDFIC 
