
#PMPNN bias

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input)
            folder="$2"
            shift # Shift past the argument value
            ;;
        --aa)
            aa_list="$2"
            shift
            ;;
        --bias)
            bias="$2"
            shift 
            ;;
        --max)
            max="$2"
            shift
            ;;
        *)
        esac
        shift
done

output_silent="input.silent"
output_bias="bias.jsonl"
mkdir -p ./slurm_logs

if ! [ -f $output_silent ]; then 
    job1=$(/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbs "$folder/*pdb" > "$output_silent")
    wait
fi
echo "Biased aminoacids: ${aa_list}"
echo "Bias introduced: ${bias}"
python3 /home/cchacon/Carlos_scripts/pMPNN/make_bias_AA.py --AA_list "$aa_list" --bias_list "$bias" --output_path "$output_bias"
wait
echo "bias json created as $output_bias"

i=1

while [[ $i -le $max ]];do

    mkdir -p "run_$i" 

    cp $output_silent run_${i}/$output_silentjid4=$(sbatch submit_af2_interfaces.sh --input_silent "$input_pmpnn")
    input_pmpnn=run_${i}/$output_silent

    echo "running pmpnn"
    job3=$(sbatch /home/cchacon/protein_design/submit_pMPNN_bias.sh --input_silent "$input_pmpnn" --n_seqs 1 --relax_cycles 0 --bias "$output_bias")

    #jid4=$(sbatch submit_af2_interfaces.sh --input_silent "$input_pmpnn")
    i=$((i+1))

done