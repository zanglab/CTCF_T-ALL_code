#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCOLN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCOLN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCRS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCTP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCTP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCTS1 
