
max=10
counter=0
for folder in ZINC* do 
    while [ $counter -le $max ]; do
        if [ ! -f $folder/$folder.pdb ]; then  
            sbatch /apps/scripts/AF2_tools/rfaa_submit_script $folder
            counter=$((counter+1))
        fi
    done 
    if [ $counter -gt $max ] ; then 
        while [ ! -f $previous_folder/$previous_folder.pdb ]; then
            sleep 120
            echo Waiting for $folder to finish
        done 
        counter=0
    fi
    previous_folder=$folder
done


