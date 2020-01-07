#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MPC1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MPC1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPDU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPDZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPEG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPG 
