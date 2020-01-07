#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LINC01556.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC01556 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINGO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINGO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINGO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINGO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPA 
