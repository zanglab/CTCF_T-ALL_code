#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MOSPD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MOSPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOV10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOV10L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOXD2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPC1 
