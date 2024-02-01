
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sc', '-sc', help='score file path', required=True)
parser.add_argument('--previous', '-pr', help='previous score file path', required=True)
args = parser.parse_args()

def check_score(actual, previous):

    actual_sc=pd.read_csv(actual,sep='\s+')
    previous_sc=pd.read_csv(previous, sep='\s+')
    if actual_sc['pae_interaction']> previous_sc['plddt_binder'] and actual_sc['plddt_binder']>threshold:
        return True 
    else:
        return False

check_score(args.sc, args.previous)

