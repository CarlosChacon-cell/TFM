for folder in *loop; do
    cd $folder
    gmx editconf -f conf.gro -d 2 -bt cubic -o box.gro
    gmx solvate -cp box.gro -cs -o solvate.gro -p topol.top
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
    #Create NVT equilibration files (all seem good)
    gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro 
    #run NVT equilibration
    gmx mdrun -v -deffnm nvt
    wait
    #Create NPT equilibration file
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    wait 
    #run NPT equilibration 
    gmx mdrun -v -deffnm npt
    wait
    #create run files for a 10ns run
    gmx grompp -f run.mdp  -c npt.gro -t npt.cpt -p topol.top -o run50ns.tpr
    wait 
    sbatch gromacs_submit_cluster_tests.sh "run50ns.tpr" "0"
    cd ..
done


