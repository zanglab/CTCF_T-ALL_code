#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o out_py2.log

time python 2_union_binding_view2M_bothside_interaction.py

