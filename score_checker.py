
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sc', '-sc', help='score file path', required=True)
parser.add_argument('--previous', '-pr', help='previous score file path', required=True)
args = parser.parse_args()

threshold=40

def check_score(actual, previous):

    actual_sc=pd.read_csv(actual,sep='\s+')
    previous_sc=pd.read_csv(previous, sep='\s+')
    if actual_sc['pae_interaction'][0]< previous_sc['pae_interaction'][0] and actual_sc['plddt_binder'][0]>threshold:
        return print(1)
    else:
        return print(0)

check_score(args.sc, args.previous)

