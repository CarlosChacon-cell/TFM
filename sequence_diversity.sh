

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input)
            input1="$2"
            shift # Shift past the argument value
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


#Setting variables to 0
i=0
best_i=0

#Copying the target to the folder 
input="${input1}.pdb"

jobx=$(pymol -c /data/cchacon/carlos/scripts/Carlos_scripts/pymol_distance_np.py --protein "$input" --peptide "$peptide" --chains "$chains" --csv "interacting_${i}" )
echo "interacting.csv created"
job1=$(python3 /data/cchacon/carlos/scripts/Carlos_scripts/fixed_trial.py --pdbs "$input" --csv interacting_${i})
echo " Residues fixed at positions ${indices}"

while [[ $i -le $max ]];do
    
    mkdir -p "./run_${i}"
    
    actual_folder="run_${i}"
    
    output_silent="$actual_folder/TFEB_run_${i}.silent"

    input_pmpnn="$actual_folder/TFEB_run_${i}.silent"
    
    output_pmpnn="$actual_folder/TFEB_run_${i}_out.silent"

    input_af2="$actual_folder/TFEB_run_${i}.silent"
        
    output_af2="${actual_folder}/TFEB_run_${i}_out_af2"

    #pdbs to silent 

    job2=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs "$input" > "$output_silent")
    echo "PDB transformed to silent with filename $actual_folder/TFEB_run_${i}.silent"

    #pmpnn_fr
    echo "starting pmpnn"
    job3=$(sbatch submit_pMPNN.sh --input_silent "$input_pmpnn" --n_seqs 1 --relax_cycles 0)
    jid3dep=`echo $job3 | awk '{print $4}'`
    echo "Submitted AF2 with jobid: $jid4dep"


    # Check if the file exists
    while [[ ! -e "$output_pmpnn" ]]; do
        echo "Waiting for the file to be created: $output_pmpnn"
        sleep 60  # Adjust the sleep duration as needed
    done

    echo "PMPNN-FR finished"

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

    i=$((i+1))

done

