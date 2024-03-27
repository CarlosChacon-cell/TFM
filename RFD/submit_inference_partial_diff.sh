#!/bin/bash

# Number of desired cpus:
#SBATCH --nodes=1
#SBATCH -p RFD
#SBATCH --open-mode=append
#SBATCH --gres=gpu:1
#SBATCH --exclusive=user
#SBATCH --cpus-per-gpu=12
#SBATCH -o slurm_logs/%j.out
#SBATCH -e slurm_logs/%j.err

export HYDRA_FULL_ERROR=1

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --output_prefix)
            output_prefix="$2"
            shift # Shift past the argument value
            ;;
        --input_pdb)
            input_pdb="$2"
            shift
            ;;
        --contigmap_descriptor)
            contigmap_descriptor="$2"
            shift
            ;;
        --hotspots_descriptor)
            hotspots_descriptor="$2"
            shift
            ;;
        --designs_n)
            designs_n="$2"
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

# Display the parsed values
machine=`hostname`
echo "Current machine $machine"diffuser.partial_T=
echo "You chose the following input values:"
echo "  - output prefix: $output_prefix"
echo "  - input PDB: $input_pdb"
echo "  - contigmap descriptor: $contigmap_descriptor"
echo "  - number of designs: $designs_n"
echo "  - number of noise steps: $noise_steps"
echo "  - noise scale: $noise_scale"

# Activate environment
source /apps/profile.d/load_all.sh
conda activate SE3nv
wait
# Run!
/apps/rosetta/RFDifussion/scripts/run_inference.py inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor"  inference.num_designs="$designs_n" diffuser.partial_T="$noise_steps" denoiser.noise_scale_ca="$noise_scale" denoiser.noise_scale_frame="$noise_scale"

echo "done"