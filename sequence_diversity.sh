

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
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

counter=0
for file in *.pdb; do
#Setting variables to 0
    i=0
    best_i=0
    
    #make a folder for each pdb file
    mkdir "design_${counter}"
    #Getting the interacting resiudes
    jobx=$(pymol -c /data/cchacon/carlos/scripts/Carlos_scripts/pymol_distance_np.py --protein "$file" --peptide "$peptide" --chains "$chains" --csv "interacting_${counter}" )
    echo "interacting.csv created"

    #while loop to generate diversity
    while [[ $i -le $max ]];do
        
        mkdir "design_${counter}/run_${i}"

        actual_folder="design_${counter}/run_${i}"
        
        input_fixed="$actual_folder/protein_${counter}.pdb"

        output_silent="$actual_folder/TFEB_run_${i}.silent"

        input_pmpnn="$actual_folder/TFEB_run_${i}.silent"
        
        output_pmpnn="$actual_folder/TFEB_run_${i}_out.silent"

        input_af2="$actual_folder/TFEB_run_${i}.silent"
            
        output_af2="${actual_folder}/TFEB_run_${i}_out_af2"

        #We save each file into the run folder 
        joby=$(pymol -c /data/cchacon/carlos/scripts/Carlos_scripts/structure_save.py --protein "$file" --i "$counter" --folder "$actual_folder")

        #we fix residues 
        job1=$(python3 /data/cchacon/carlos/scripts/Carlos_scripts/fixed_trial.py --pdbs "$input_fixed" --csv interacting_${counter})
        echo " Residues fixed at positions ${indices}"

        #pdbs to silent 

        job2=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs "$input_fixed" > "$output_silent")
        echo "PDB transformed to silent with filename $actual_folder/TFEB_run_${i}.silent"

        #pmpnn_fr
        echo "starting pmpnn"
        job3=$(sbatch submit_pMPNN.sh --input_silent "$input_pmpnn" --n_seqs 1 --relax_cycles 1)
        jid3dep=`echo $job3 | awk '{print $4}'`
        echo "Submitted pMPNN with jobid: $jid3dep"


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
    counter=$((counter+1))
done
