#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOR1AIP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR1AIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOX 
