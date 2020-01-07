#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLCXD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCXD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCXD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLD2 
