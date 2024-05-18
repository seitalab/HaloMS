import json
import matplotlib.pyplot as plt
import os
import py3Dmol

import statistics
import math

import pandas as pd
import glob
        
        
def get_score(inputdirectorypath, outputdirectorypath , seed):
    input_files = sorted(glob.glob(inputdirectorypath + '/*'))
    project_name = inputdirectorypath.split('/')[-1]
    print(project_name)
    output_score = []
    output_score_csv = 'output/' + project_name + '_all_score.csv'
    
    for file in input_files:    
        file_path = outputdirectorypath + '/' + outputdirectorypath.split('/')[-1]
        ranking_debug_path = '{}_unrelaxed_rank_1_model_1_scores.json'.format(file_path)
        print(ranking_debug_path)
        
        tmp_path = '{}_template_domain_names.json'.format(file_path)
        pro_a = file_path.split('__')[-3]
        pro_b = file_path.split('__')[-2]
        print(pro_a)
        
        a_tmp = []
        b_tmp = []
        if os.path.isfile(tmp_path):
            with open(tmp_path) as f:
                jsn = json.load(f)
                # print('A template :', jsn['A'])
                # print('B template :', jsn['B'])
                a_tmp = jsn['A']
                b_tmp = jsn['B']


        if os.path.isfile(ranking_debug_path):
            with open(ranking_debug_path) as f:
                jsn = json.load(f)
            output_score.append([project_name, pro_a, pro_b,
                                         jsn['iptm']*0.8+jsn['ptm']*0.2, jsn['iptm'], jsn['ptm'], a_tmp, b_tmp])

        else:
            output_score.append([project_name, pro_a, pro_b, 'error', 'error', 'error', 'error', 'error'])

    df = pd.DataFrame(output_score, columns=['project_name', 'Protein_A', 'Protein_B', 'model_1_multimer_score(iptm+ptm)', 'model_1_multimer_score(iptm)', 'model_1_multimer_score(ptm)', 'A template', 'B template'])
    
    df.to_csv(output_score_csv)
    print('output score : ' + output_score_csv)


# Please enter the your data path in XXXX, YYYY
INPUT_FOLDER='XXXX'
OUTPUT_FOLDER='YYYY'
SEED=0
get_score(INPUT_FOLDER, OUTPUT_FOLDER, SEED)