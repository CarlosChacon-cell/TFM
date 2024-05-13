
for protein in *pdb ; do 

    protein_name=$(echo "$protein" | sed "s/.pdb//")
    cd "$protein_name"
    if [[ ! -f "run2ns.tpr" ]]; then
        gmx editconf -f conf.gro -d 2 -bt cubic -o box.gro
        gmx solvate -cp box.gro -cs -o solvate.gro -p topol.top 
        gmx grompp -f ions.mdp -c solvate.gro -p topol.top -o ions.tpr
        echo '13' | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
        gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr
        gmx mdrun -v -deffnm em
        wait
        gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
        gmx mdrun -v -deffnm nvt
        wait
        gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
        gmx mdrun -v -deffnm npt
        wait
        gmx grompp -f run.mdp  -c npt.gro -t npt.cpt -p topol.top -o run2ns.tpr
        
        for i in {1..5}; do 
            mkdir -p run_${i}
            cp run2ns.tpr run_$i/.
            cp *itp run_$i/.
            cp *top run_$i/.
            cd run_$i
            sbatch gromacs_submit_cluster.sh "run2ns.tpr" "0"
            cd ..
        done
    wait
        # Wait until a .gro file is created
        while [ ! -f "run_5/production_run.part0001.gro" ]; do
            sleep 10  # Wait for 10 seconds before checking again
        done
        for i in {1..5}; do 
            cd run_${i}
            echo "1 0" | gmx trjconv -s run2ns.tpr -f production_run.part0001.xtc -o noPBC.xtc -pbc nojump -center
            cd ..
        done
    fi
    cd ..
done 

