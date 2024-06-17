'''
RFAA needs as input a yaml file. This codes generates that yaml file.
-->Inputs:
    --fasta: Path to the fasta file of the target protein. It can be one or two fastas, separated by a whitespace
    --comp: Path to the compound file. This compound must be in sdf format, and right now only a compound is admitted
    --folder: Path to the folder in which the yaml file is going to be written. The yaml structure is prepared to store the yaml file in the folder
      from which RFAA is going to be run, with a structure like: 
      
      Parent_directory/fastas
      Parent_directory/compounds
      Parent_directory/running_folder

-->Output:
    input.yaml. A yaml folder with the needed structure to run RFAA

'''
import os 
import argparse


def get_file_type(file_path):
    _, file_extension = os.path.splitext(file_path)
    file_extension=file_extension.split('.')[1]
    return file_extension.lower()

def get_file_name(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

def create_yaml_content(fasta_files, comp_path, comp_name):
    # Define the content structure
    if len(fasta_files)==2:
        defaults = [
            {
                'base_cluster': None,
                'job_name': comp_name,
                'protein_inputs': {
                    'A': {
                        'fasta_file': fasta_files[0]
                    },
                    'B': {
                        'fasta_file': fasta_files[1]
                    }
                },
                'sm_inputs': {
                    'C': {
                        'input': comp_path,
                        'input_type': "sdf"  # Remove single quotes, just in case
                    }
                }
            }
        ]
    elif len(fasta_files)==1:
        defaults = [
            {
                'base_cluster': '',
                'job_name': comp_name,
                'protein_inputs': {
                    'A': {
                        'fasta_file': fasta_files[0]
                    }
                },
                'sm_inputs': {
                    'C': {
                        'input': comp_path,
                        'input_type': "sdf"  # Remove single quotes, just in case
                    }
                }
            }
        ]

    return {'defaults': defaults}

def write_yaml_content(yaml_content, filename, numb_args):
    # Create YAML formatted string manually
    yaml_string = "defaults:\n"
    if numb_args==2:
        for item in yaml_content['defaults']:
            yaml_string += "- base_cluster\n"
            yaml_string += 'job_name: "{}"\n'.format(item['job_name'])
            yaml_string += "protein_inputs:\n"
            yaml_string += "      A:\n"
            yaml_string += "        fasta_file: ../{}\n".format(item['protein_inputs']['A']['fasta_file'])
            yaml_string += "      B:\n"
            yaml_string += "        fasta_file: ../{}\n".format(item['protein_inputs']['B']['fasta_file'])
            yaml_string += "sm_inputs:\n"
            yaml_string += "      C:\n"
            yaml_string += "        input: ../{}\n".format(item['sm_inputs']['C']['input'])
            yaml_string += '        input_type: \"{}\"\n'.format(item['sm_inputs']['C']['input_type'])
    elif numb_args==1:
          for item in yaml_content['defaults']:
            yaml_string += "- base_cluster\n"
            yaml_string += 'job_name: "{}"\n'.format(item['job_name'])
            yaml_string += "protein_inputs:\n"
            yaml_string += "      A:\n"
            yaml_string += "        fasta_file: ../{}\n".format(item['protein_inputs']['A']['fasta_file'])
            yaml_string += "sm_inputs:\n"
            yaml_string += "      B:\n"
            yaml_string += "        input: ../{}\n".format(item['sm_inputs']['C']['input'])
            yaml_string += '        input_type: \"{}\"\n'.format(item['sm_inputs']['C']['input_type'])      
    # Write the manually created YAML string to the file
    with open(filename, 'w') as yaml_file:
        yaml_file.write(yaml_string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', help='protein fasta file path', required=False, nargs='+')
    parser.add_argument('--comp', help='comp file path', required=True)
    parser.add_argument('--folder', help='folder in which you want to save the yaml_file ', required=True)
    args = parser.parse_args()

    comp_path = args.comp
    comptype = get_file_type(comp_path)
    comp_name = get_file_name(comp_path)
    numb_args=len(args.fasta)
    if args.fasta:
        filename = f'{args.folder}/input.yaml'
        yaml_content = create_yaml_content(args.fasta, args.comp, comp_name)
        write_yaml_content(yaml_content, filename,numb_args)
