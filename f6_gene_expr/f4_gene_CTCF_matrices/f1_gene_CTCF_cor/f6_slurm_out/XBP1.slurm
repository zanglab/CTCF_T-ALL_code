#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o XBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i XBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i XG 
