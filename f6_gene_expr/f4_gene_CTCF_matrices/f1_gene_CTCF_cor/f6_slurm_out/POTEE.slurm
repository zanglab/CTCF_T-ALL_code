#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POTEE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POTEJ 
