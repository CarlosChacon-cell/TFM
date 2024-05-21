
#This is just to build the ymla files neccesary for the RFAA from the PDB/fasta and compfile
import os 
import argparse
import re 
import yaml


def pdb_to_fasta(pdb_file):
    sequence = ""

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[13:15] == "CA":
                # Extract the amino acid code
                aa = line[17:20].strip()
                sequence+=aa

                # Append the amino acid to the sequence
                        
    return sequence

# def create_yaml_content(fasta_files, comp_file,comp_name):
#     # Generate protein_inputs dictionary based on fasta files
#     protein_inputs = {}
#     if len(fasta_files)==2:
#         protein_inputs= {
#             'A': {
#                 'fasta_file':f'../{fasta_files[0]}'
#             },
#             'B': {
#                 'fasta_file':f'../{fasta_files[1]}'
#             }
#         }
#     else:
#         protein_inputs= {
#            'A': {
#                 'fasta_file':f'../{fasta_files[0]}'
#             }                      
#             }
    
#     # Hardcoded sm_inputs dictionary

#     sm_inputs = {
#         'C': {
#             'input': f'../{comp_file}',
#             'input_type': "\"sdf\""
#         }
#     }
    
#     # YAML structure
#     yaml_content = {
#         'defaults': ['base_cluster'],
#         'job_name': f'"{comp_name}"',
#         'protein_inputs': protein_inputs,
#         'sm_inputs': sm_inputs
#     }
#     return yaml_content

def get_file_type(file_path):
    _, file_extension = os.path.splitext(file_path)
    file_extension=file_extension.split('.')[1]
    return file_extension.lower()

def get_file_name(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]




# if args.fasta:
#     filename=(f'{args.folder}/input.yaml')
#     yaml_content = create_yaml_content(args.fasta, comp_path, comp_name)
#     with open(filename, 'w') as yaml_file:
#         yaml.dump(yaml_content, yaml_file, default_flow_style=False)

# # Read YAML data
# with open(filename, 'r') as file:
#     data = yaml.safe_load(file)

# print(data['job_name'])
# # Modify fields

# data['job_name'] = data["job_name"].strip("'")
# data['sm_inputs']['C']['input_type']=data['sm_inputs']['C']['input_type'].strip('')

# # Write modified YAML data
# with open(filename, 'w') as file:
#     yaml.dump(data, file)


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
    parser.add_argument('--pdb', help='protein pdb file path', required=False)
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
