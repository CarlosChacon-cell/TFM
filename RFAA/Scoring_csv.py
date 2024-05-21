import argparse
import pandas as pd
from Scoring import extract_scoring
import glob 

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--folder', '-f')
    args=parser.parse_args()

    filename=glob.glob(f'{args.folder}/*_aux.pt')[0]
    mean_plddt, pae_interaction, pae=extract_scoring(filename)
    score_dict={
        'Comp_name':args.folder,
        'mean_plddt': [round(mean_plddt,3)],
        'pae_interaction': [round(pae_interaction,3)]
    }
    score_df=pd.DataFrame(score_dict)
    try:
        old_df=pd.read_csv('LigandSearchMetrics.csv', index_col=0)
        merge_df=pd.concat([old_df, score_df])
        merge_df.to_csv('LigandSearchMetrics.csv')
    except:
        score_df.to_csv('LigandSearchMetrics.csv')
        
