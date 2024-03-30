
:'This is the key code for the (sadly not working) self-updating pMPNN cycle
Future works must be done, specially in one direction as point out by a Baker paper:
Use the sequence conservation of a protein to improve its fitness. Keeping the 50% residues they
could improve the function of TEV protease. From now on, useful comments'


:'This code allows you to change the sequence of a binder in order to improve its binding efficiency
You can keep (and its highly recommended) the residues important for the interaction'

:'I think the sequence diversity code is more promising so better check that one (pretty similar) before this'

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input)
            input1="$2" #The name of the protein and the binder without the pdb
            shift # Shift past the argument value
            ;;
        --indices) 
            indices="$2" #Residues to be kept 
            shift 
            ;;
        --max)
            max="$2" #Number of cycles that will be done
            shift 
            ;;
        --peptide)
            peptide="$2" #Peptide chain
            shift
            ;;
        --chains)
            chains="$2" #Target chain (all in one chain preferably)
            shift
            ;;
        --name)
            name="$2" #Name of the run
            shift 
            ;;
        *)
        esac
        shift
done

i=0
best_i=0

#input1 is the input without the pdb termination so it can be handle easily



while [[ $i -le $max ]];do
    
    #We name all needed variables 

    mkdir -p "./run_${i}"
    
    input="$input1.pdb"

    actual_folder="run_${i}"

    cp "$input" "$actual_folder/"
    
    input_fixed="$actual_folder/$input"
    
    output_silent="$actual_folder/${name}_run_${i}.silent"

    input_pmpnn="$actual_folder/${name}_run_${i}.silent"
    
    output_pmpnn="$actual_folder/${name}_run_${i}_out.silent"

    input_af2="$actual_folder/${name}_run_${i}.silent"
        
    output_af2="${actual_folder}/${name}_run_${i}_out_af2"

    #fixing residues
    job1=$(python3 /data/carlos/scripts/Carlos_scripts/fixed_trial.py --pdbs "$input_fixed" --indices "$indices" --csv interacting_${best_i})
    echo " Residues fixed at positions ${indices}"
    
    #pdbs to silent 

    job2=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs "$input_fixed" > "$output_silent")
    echo "PDB transformed to silent with filename $actual_folder/${name}_run_${i}.silent"
   

    #pmpnn_fr
    echo "starting pmpnn"
    job3=$(sbatch submit_pMPNN.sh --input_silent "$input_pmpnn" --n_seqs 1 --relax_cycles 1)

    # Check if the file exists
    while [[ ! -e "$output_pmpnn" ]]; do
        echo "Waiting for the file to be created: $output_pmpnn"
        sleep 60  # Adjust the sleep duration as needed
    done

    echo "PMPNN-FR finished, outfile is $actual_folder/${name}_run_${i}_pmpnn_out.silent"
    #submit af2

    job4=$(sbatch submit_af2_interfaces.sh --input_silent "$input_af2")
    jid4dep=`echo $job4 | awk '{print $4}'`
    echo "Submitted AF2 with jobid: $jid4dep"

    
    # Specify the file to wait for
        AF2sc="${output_af2}.sc"
        # Check if the file exists
        while [[ ! -e "$AF2sc" ]]; do
            echo "Waiting for the file to be created: $AF2sc"
            sleep 60  # Adjust the sleep duration as needed
        done

        # Continue with the rest of your code after the file is created
        echo "AF2 finished, continuing with the script"

    if [ "$i" -ne 0 ];then #If we are no in the first cycle we check the pae..
        #Check if the new generated protein is better than the previous best
        result=$(python3 /data/carlos/scripts/Carlos_scripts/score_checker.py --sc "run_${i}/${name}_run_${i}_out_af2.sc" --previous "run_${best_i}/${name}_run_${best_i}_out_af2.sc")
        if [ "$result" -eq 1 ] ; then #If it is better
            #extract the pdb 
            job5=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentextract "${output_af2}.silent")
            echo "Silent extracted as ${input1}_dldesign_0_af2pred.pdb"
            best_i=$i
            #Get the distances 
            job6=$(pymol -c /data/carlos/scripts/Carlos_scripts/pymoltrial.py --protein "${input1}_dldesign_0_af2pred.pdb" --peptide "$peptide" --chains "$chains" --csv "/interacting_${i}" --i "$i")
            echo "Distances computed" 
            input1="protein_${i}"
            echo "new best score is cycle ${best_i}"
            #Save the structure
            job7=$(pymol -c /data/carlos/scripts/Carlos_scripts/structure_save.py --protein "${input1}_dldesign_0_af2pred.pdb" --i "$i")
            input1="protein_${i}"
            echo "best score is cycle ${best_i}"

        else #If it wotse, just save teh structure but NOT update the input for the next cycle
            job8=$(pymol -c /data/carlos/scripts/Carlos_scripts/structure_save.py --protein "${input1}.pdb" --i "$i")
            input1="protein_${best_i}"
            echo "best score is cycle ${best_i}"
        fi
        
    else #If we are in the first cycle, just save the structure 

        job9=$(pymol -c /data/carlos/scripts/Carlos_scripts/structure_save.py --protein "${input1}.pdb" --i "$i")
        input1="protein_${i}"
        echo "best score is cycle ${best_i}"

    fi 

    i=$((i+1))

done

