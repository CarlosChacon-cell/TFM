

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input)
            input1="$2"
            shift # Shift past the argument value
            ;;
        --indices)
            indices="$2"
            shift 
            ;;
        --max)
            max="$2"
            shift 
            ;;
        --peptide)
            peptide="$2"
            shift
            ;;
        --chains)
            chains="$2"
            shift
            ;;
        *)
        esac
        shift
done

i=0
best_i=1

#input1 is the input without the pdb termination so it can be handle easily

input="$input1.pdb"

#pattern of the fixed
# pattern="^REMAR*"
# Specify the input file

while [[ $i -le $max ]];do
    
    mkdir -p "./run_${i}"
    
    input="$input1.pdb"

    actual_folder="run_${i}"

    cp "$input" "$actual_folder/"
    
    input_fixed="$actual_folder/$input"
    
    output_silent="$actual_folder/TFEB_run_${i}.silent"

    input_pmpnn="$actual_folder/TFEB_run_${i}.silent"
    
    output_pmpnn="$actual_folder/TFEB_run_${i}_pmpnn_out.silent"

    input_af2="$actual_folder/TFEB_run_${i}_pmpnn.silent"
        
    output_af2="${actual_folder}/TFEB_run_${i}_pmpnn_out_af2"

    #fixing residues
    job1=$(python3 /data/carlos/scripts/fixed_trial.py --pdbs "$input_fixed" --indices "$indices" --csv ./distance_outputs/interacting_${best_i})
    echo " Residues fixed at positions ${indices}"
    
    #pdbs to silent 

    job2=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs "$input_fixed" > "$output_silent")
    echo "PDB transformed to silent with filename $actual_folder/TFEB_run_${i}.silent"
   

    #pmpnn_fr
    echo "starting pmpnn"
    job3=$(python3 -u /apps/rosetta/dl_binder_design/mpnn_fr/dl_interface_design.py -silent "$input_pmpnn" -checkpoint_path "/apps/rosetta/dl_binder_design/mpnn_fr/ProteinMPNN/vanilla_model_weights/v_48_030.pt" -outsilent "$output_pmpnn")
    echo "PMPNN-FR finished, outfile is $actual_folder/TFEB_run_${i}_pmpnn_out.silent"

    #pmpnn_fr becomes silly if the pdb from which you want to create a new sequence has already been inside pmpnn.
    #The next line is thought to stop the cycle if the new file is not created

    if [ ! -e "$output_pmpnn" ]; then
        echo "Error: The required file '$actual_folder/TFEB_run_${i}_pmpnn_out.silent' does not exist. Exiting."
        exit 1  # Exit with a non-zero status to indicate an error
    fi

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

    #extract the pdb 

    job5=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentextract "${output_af2}.silent")
    echo "Silent extracted as ${input1}_dldesign_0_cycle1_af2pred.pdb"

    #now we check the polar interactions to fix it 
    job6=$(pymol /data/carlos/scripts/pymoltrial.py --protein "${input1}_dldesign_0_cycle1_af2pred.pdb" --peptide "$peptide" --chains "$chains" --csv "./distance_outputs/interacting_${i}")
    jid6dep=`echo $job6 | awk '{print $4}'`
    echo "Distances computed"

    mv "${input1}_dldesign_0_cycle1_af2pred.pdb" "output_${i}.pdb"
    if [ "$i" -gt 0 ];then
        if $(python3 /data/carlos/scripts/score_checker.py --sc "TFEB_run_${i}_pmpnn_out_af2.sc" --previous "TFEB_run_${best_i}_pmpnn_out_af2.sc"); then
             best_i=$i
             input1="output_${i}"
             echo "new best score is cycle ${i}"
        fi
    fi 

    i=$((i+1))

done

