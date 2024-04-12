#! /bin/bash
source /apps/profile.d/load_all.sh

# Master script for deep searches in microruns

# 1: Get all info needed
## Parse command-line arguments

#Set defaults

partial_diff="False"
noise_steps=20
pmp_nseqs=1
rfd_ndesigns=8
pmp_relax_cycles=1
noise_scale=1
#contigs_getter
rfd_contigs=$(contigs_getter_pd.py --file $input)



while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input)
            input="$2"
            shift # Shift past the argument value
            ;;
        --template)
            template="$2"
            shift 
            ;;
        --max_threads)
            max="$2"
            shift 
            ;;    
        --rfd_contigs) #   Automatic for partial diffusion
            rfd_contigs="$2"
            shift
            ;;    
        --rfd_hotspots)
            rfd_hotspots="$2"
            shift 
            ;;    
        --rfd_ndesigns)
            rfd_ndesigns="$2"
            shift
            ;;    
        --pmp_nseqs)
            pmp_nseqs="$2"
            shift
            ;;    
        --pmp_relax_cycles)
            pmp_relax_cycles="$2"
            shift
            ;;   
        --partial_diff)
            partial_diff="$2"
            shift 
            ;; 
        --noise_steps)
            noise_steps="$2"
            shift
            ;;
        --noise_scale)
            noise_scale="$2"
            shift
            ;;
        
        *)
        echo "Unknown option: $1"
            exit 1
            ;;


    
    esac
    shift # Shift past the current argument
done

# Get ready

mkdir -p ./slurm_logs
mkdir -p ./output
mkdir -p ./jsons

last_run_folder=$(ls -d "./output/run_"* 2>/dev/null | sort -V | tail -n 1 | sed 's#./output/run_##')

if [[ -n "$last_run_folder" ]]; then
    i="$last_run_folder"
else
    i=0
fi




# RUN
while true; do
    i=$((i+1))
    echo "Starting cycle $i"
    mkdir -p ./output/run_$i

    # Set a default value for previous if i is smaller than max (initial cycles)
    if [ "$i" -le $max ]; then
        previous=0
    else
        previous=$((i - $max))
    fi
    
    # Set input, output and helping filenames
    output_rfd="output/run_$i/run_${i}_design"
    input_pymol="output/run_$i/run_${i}_design_0.pdb"
    output_pymol="output/run_$i/run_${i}_template_aligned.pdb"
    input_sub="output/run_$i/"
    output_sub="output/run_$i/"
    output_silent="run_${i}_input.silent"
    input_pmp="output/run_$i/run_${i}_input.silent"
    waitfor="output/run_${previous}/run_${previous}_done"
    end_touch="output/run_$i/run_${i}_done"

    # Skip waiting when previous is 0
    if [ "$previous" -ne 0 ]; then
        # Wait until the condition is met
        while [ ! -e "$waitfor" ]; do
            echo "Waiting for previous jobs to complete: $waitfor"
            sleep 60
        done
    # Now that the condition is met or previous is 0, proceed to the following code
    fi
    
    # 2: Run RFD without dependencies and with partial diffusion
    #Note on how it works: It seems you just have to add 
    #'contigmap.contigs=[len-len/0 B1-len]' (I think is 0 index so len-1)
    #You specify that the chain A is the one to modify and that the chain B must be keep intact
    #You have to also add the diffuser.partial_T=20 (recommended, 2/5 of the total numbers of step for noising )
    if [ $partial_diff = "True" ]; then
        jid1=$(sbatch submit_inference_partial_diff.sh --output_prefix "$output_rfd" --input_pdb "$input" --contigmap_descriptor "$rfd_contigs"  --designs_n "$rfd_ndesigns" --noise_steps "$noise_steps" --noise_scale "$noise_scale")
        jid1dep=`echo $jid1 | awk '{print $4}'`
        echo "Submitted RFD with jobid: $jid1dep"
    else 
        jid1=$(sbatch submit_run_inference.sh --output_prefix "$output_rfd" --input_pdb "$input" --contigmap_descriptor "$rfd_contigs" --hotspots_descriptor "$rfd_hotspots"  --designs_n "$rfd_ndesigns")
        jid1dep=`echo $jid1 | awk '{print $4}'`
        echo "Submitted RFD with jobid: $jid1dep"
    fi

    #Maybe add that if there is no partial diff go to the RFD route ??
    # 3: Run Pymol align + chain substitute + silent creation splitting into 4 chunks
    jid2=$(sbatch --dependency=afterany:$jid1dep submit_pymolalign_subtitutechain_silentgroup.sh --pymol_template "$template" --pymol_reference "$input_pymol" --pymol_output "$output_pymol" --sub_input "$input_sub" --sub_output "$output_sub" --silent_output "$output_silent" --rfd_ndesigns "$rfd_ndesigns")
    jid2dep=`echo $jid2 | awk '{print $4}'`
    echo "Submitted align-substitute-silent with jobid $jid2dep waiting for $jid1dep to finish"

    # 4: Run pMPNN
    jid3=$(sbatch --dependency=afterany:$jid2dep submit_pMPNN.sh --input_silent "$input_pmp" --n_seqs "$pmp_nseqs" --relax_cycles "$pmp_relax_cycles")
    jid3dep=`echo $jid3 | awk '{print $4}'`
    echo "Submitted pMPNN for $split_silent with jobid: $jid3dep waiting for $jid2dep to finish"

    # 5: Run AF2
    jid4=$(sbatch --dependency=afterany:$jid3dep submit_af2_interfaces.sh --input_silent "$input_pmp")
    jid4dep=`echo $jid4 | awk '{print $4}'`
    echo "Submitted AF2 for $split_silent with jobid: $jid4dep waiting for $jid3dep to finish"

    # 6: All done, mark the end of the microrun for next rounds of the while loop to be allowed
    jid5=$(sbatch --dependency=afterany:$jid4dep submit_finish_microrun.sh --output "$end_touch")
    jid5dep=`echo $jid5 | awk '{print $4}'`
    echo "Submitted final touch fo run $i with jobid: $jid5dep waiting for $jid4dep to finish"

    #7: Move all json to a json folder
    wait 
    mv *json jsons/.

done
