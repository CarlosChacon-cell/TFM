while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input)
            input1="$2"
            shift # Shift past the argument value
            ;;
        *)
        esac
        shift
done

i=0

while [[ $i -le 10 ]];do

    mkdir "run_$i"

    cp "design_1034_out.silent" ./"run_$i"/

    cd ./"run_$i"

    mkdir "slurm_logs"

    job4=$(sbatch submit_af2_interfaces.sh --input_silent "$input1")
    jid4dep=`echo $job4 | awk '{print $4}'`
    echo "Submitted AF2 with jobid: $jid4dep"

    AF2sc="run_${i}/design_1034_out_af2.sc"
    # Check if the file exists
    while [[ ! -e "$AF2sc" ]]; do
        echo "Waiting for the file to be created: $AF2sc"
        sleep 60  # Adjust the sleep duration as needed
    done
    
    cd ..

    i=$((i+1))

done