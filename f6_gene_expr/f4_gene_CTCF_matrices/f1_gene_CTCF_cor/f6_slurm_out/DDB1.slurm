#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDHD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDI1 
