#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DIAPH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DIAPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIAPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIAPH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DICER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIDO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIEXF 
