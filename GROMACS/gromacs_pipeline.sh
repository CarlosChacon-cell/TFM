

#GROMACS PIPELINE

for i in {10..44}; do
    cd design_$i
    file="*.pdb"
    #Transform to gmx
    gmx pdb2gmx -f $file -ignh -ff "charmm36" -water "tip3p" #Ignore hydrogens, only important with protonation states

    if [ ! -e "conf.gro" ];then
        cd ..
        continue
    fi
    #Change the size of the box, 2 A at each side 
    gmx editconf -f conf.gro -bt triclinic -d 1.0 -o box.gro
    #Add waters to the box
    gmx solvate -cp box.gro -cs -o solvate.gro -p topol.top
    #Create the ions file based on a mdp (I think any .mdp is valid)
    gmx grompp -f ions.mdp -c solvate.gro -p topol.top -o ions.tpr
    #Add ions to the system to neutralize charges
    echo "13"|gmx genion -s ions.tpr -o solvate.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 #Conc to 0.15M to mimick in vivo conditions
    #Create minimization tpr file from the em.mdp (this time it is important) 
    gmx grompp -f em.mdp -c solvate.gro -p topol.top -o em.tpr
    #Run minimization
    gmx mdrun -v -deffnm em
    #Create NVT equilibration files (all seem good)
    gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro
    #run NVT equilibration
    gmx mdrun -v -deffnm nvt
    #Create NPT equilibration file
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    #run NPT equilibration 
    gmx mdrun -v -deffnm npt
    #create run files for a 10ns run
    gmx grompp -f run.mdp  -c npt.gro -t npt.cpt -p topol.top -o run10ns.tpr
    wait 
    #PRODUCTION RUN 20 mins approx
    sbatch gromacs_submit_cluster.sh "run10ns.tpr" "0,1,2,3"
    #gmx mdrun -deffnm run10ns_2 -nb gpu #-nb gpu makes sure you ran the simulation on the GPU
    #Remove PBC
    while [[ ! -e  "production_run.part0001.gro" ]]; do
        echo "Waiting for the file to be created"
        sleep 60  # Adjust the sleep duration as needed
    done
    echo '1 0' | gmx trjconv -s run10ns.tpr -f production_run.part0001.xtc -o production_run.part0001_NoPBC.xtc -pbc mol -center
    wait
    #get the rmsd 
    echo '4 4'| gmx rms -s run10ns.tpr -f production_run.part0001_NoPBC.xtc -o rmsd.xvg -tu ns
    wait
    #get the radius gyration
    echo '1'|gmx gyrate -s run10ns.tpr -f production_run.part0001_NoPBC.xtc -o gyrate.xvg
    wait
    cd ..
done 







        
