
import torch 
import argparse

def extract_scoring(filename):
    scoring_dict=torch.load(filename, map_location='cpu')
    mean_plddt=scoring_dict['mean_plddt']
    pae_interaction=scoring_dict['pae_inter']
    pae=scoring_dict['pae']
    return mean_plddt, pae_interaction, pae


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--file', '-f', help='File to be analysed')
    args=parser.parse_args()
    filename=args.file
    mean_plddt,pae_interaction,pae=extract_scoring(filename)
    print('The mean plddt of the prediction is: ', mean_plddt)
    print('\nThe global pae_interaction is: ', pae_interaction)
    