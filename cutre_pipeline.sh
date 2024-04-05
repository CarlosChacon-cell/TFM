
for folder in *; do 
    cd $folder
    silent="*silent"
    python3 -u /apps/rosetta/dl_binder_design/af2_initial_guess/predict_cutre.py -silent $silent -outsilent out.silent -scorefilename scores.sc  -checkpoint_name score.point
    wait
    silentextract $silent
    wait
    for protein in *pdb; do
        pymol -c /emdata/cchacon/scripts/close_residues.py --protein $protein
        wait
        python3 /emdata/cchacon/scripts/pae_v2.py --protein $protein 
    done 
    cd ..
done 