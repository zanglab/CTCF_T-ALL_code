#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GON4L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GON4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GORAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GORASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GORASP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOSR1 
