

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
best_i=0

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
    
    output_pmpnn="$actual_folder/TFEB_run_${i}_out.silent"

    input_af2="$actual_folder/TFEB_run_${i}.silent"
        
    output_af2="${actual_folder}/TFEB_run_${i}_out_af2"

    #fixing residues
    job1=$(python3 /data/carlos/scripts/Carlos_scripts/fixed_trial.py --pdbs "$input_fixed" --indices "$indices" --csv ./interacting_${best_i})
    echo " Residues fixed at positions ${indices}"
    
    #pdbs to silent 

    job2=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs "$input_fixed" > "$output_silent")
    echo "PDB transformed to silent with filename $actual_folder/TFEB_run_${i}.silent"
   

    #pmpnn_fr
    echo "starting pmpnn"
    #job3=$(sbatch  submit_pMPNN.sh --input_silent "$input_pmpnn")
    job3=$(sbatch submit_pMPNN.sh --input_silent "$input_pmpnn" --n_seqs 1 --relax_cycles 0)

    # Check if the file exists
    while [[ ! -e "$output_pmpnn" ]]; do
        echo "Waiting for the file to be created: $output_pmpnn"
        sleep 60  # Adjust the sleep duration as needed
    done

    echo "PMPNN-FR finished, outfile is $actual_folder/TFEB_run_${i}_pmpnn_out.silent"
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
    echo "Silent extracted as ${input1}_dldesign_0_af2pred.pdb"

    if [ "$i" -ne 0 ];then
        result=$(python3 /data/carlos/scripts/Carlos_scripts/score_checker.py --sc "run_${i}/TFEB_run_${i}_out_af2.sc" --previous "run_${best_i}/TFEB_run_${best_i}_out_af2.sc")
        if [ "$result" -eq 1 ] ; then
            best_i=$i
            job6=$(pymol /data/carlos/scripts/Carlos_scripts/pymoltrial.py --protein "${input1}_dldesign_0_af2pred.pdb" --peptide "$peptide" --chains "$chains" --csv "/interacting_${i}" --i "$i")
            echo "Distances computed" 
            input1="output_${i}"
            echo "new best score is cycle ${best_i}"
        else
            job7=$(pymol /data/carlos/scripts/Carlos_scripts/structure_save.py --protein "${input1}.pdb" --i "$i")
            input1="output_${i}"
            echo "best score is cycle ${best_i}"
        fi
        
    else

        job8=$(pymol /data/carlos/scripts/Carlos_scripts/structure_save.py --protein "${input1}.pdb" --i "$i")
        input1="output_${i}"
        echo "best score is cycle ${best_i}"

    fi 

    i=$((i+1))

done

