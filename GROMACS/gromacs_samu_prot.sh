
#Create the ions file based on a mdp (I think any .mdp is valid)
gmx grompp -f ions.mdp -c solvate.gro -p topol.top -o ions.tpr
#Add ions to the system to neutralize charges
echo "13"|gmx genion -s ions.tpr -o solvate.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 #Conc to 0.1M to mimick mitochondria conditions
#Create minimization tpr file from the em.mdp (this time it is important) 
gmx grompp -f em.mdp -c solvate.gro -p topol.top -o em.tpr -maxwarn 1
wait 
#Run minimization
gmx mdrun -v -deffnm em
wait 
echo '11 0' | gmx energy -f em.edr -o potential.xvg
#Create NVT equilibration files (all seem good)
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro 
#run NVT equilibration
gmx mdrun -v -deffnm nvt
wait
echo '17 0' | gmx energy -f nvt.edr -o temperature.xvg
#Create NPT equilibration file
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
wait 

#run NPT equilibration 
gmx mdrun -v -deffnm npt
wait
echo '25 0' | gmx energy -f npt.edr -o density.xvg
echo '19 0' | gmx energy -f npt.edr -o pressure.xvg
#create run files for a 10ns run
gmx grompp -f run.mdp  -c npt.gro -t npt.cpt -p topol.top -o run20ns.tpr
wait 

#prepare all files

mkdir plots 

for file in *xvg; do
    python3 /home/cchacon/cchacon/carlos/scripts/Carlos_scripts/GROMACS/remove_lines.py --file $file
    python3 /home/cchacon/cchacon/carlos/scripts/Carlos_scripts/GROMACS/gromacs_visual.py --file temperature.xvg --mode T --folder plots --i 1
    python3 /home/cchacon/cchacon/carlos/scripts/Carlos_scripts/GROMACS/gromacs_visual.py --file pressure.xvg --mode P --folder plots --i 1
    python3 /home/cchacon/cchacon/carlos/scripts/Carlos_scripts/GROMACS/gromacs_visual.py --file density.xvg --mode D --folder plots --i 1
    python3 /home/cchacon/cchacon/carlos/scripts/Carlos_scripts/GROMACS/gromacs_visual.py --file potential.xvg --mode E --folder plots --i 1
done

