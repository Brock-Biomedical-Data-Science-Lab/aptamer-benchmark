#!/bin/bash
cd "$(pwd)"

# read CUDA_VISIBLE_DEVICES settings
echo "Using GPU $CUDA_VISIBLE_DEVICES"

# === Energy Minimization ===
gmx grompp -f step4.0_minimization.mdp -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -ntmpi 5 -ntomp 8 -nb gpu

# === NVT ===
gmx grompp -f step4.1_equilibration.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu

# === Production MD (NPT) ===
gmx grompp -f step5_production.mdp -c nvt.gro -p topol.top -n index.ndx -o md.tpr
gmx mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu -v