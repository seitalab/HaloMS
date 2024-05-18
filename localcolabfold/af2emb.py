import numpy as np
import glob
import os
import pandas as pd


# Please enter the your data path in XXXX, YYYY, ZZZZ

project_name = 'XXXX'
out_dir_base = 'YYYY' + project_name + '/'
c_name = 'ZZZZ'
result_file = c_name + project_name + '__'
out_dir = glob.glob(out_dir_base + '*')


full_npy=[]
sample_id = []
for r in out_dir:
    r_id = r.replace(out_dir_base, '')
    sample_id.append(r_id.replace(c_name, ''))
    if os.path.isfile(r + '/{}_single_repr_rank_001_alphafold2_ptm_model_1_seed_000.npy'.format(r_id)):
        em = np.load(r + '/{}_single_repr_rank_001_alphafold2_ptm_model_1_seed_000.npy'.format(r_id))        
    else:
        em = [['error']*256]
    
    full_npy.append(em[0])
    
full_npy_data = pd.DataFrame(full_npy, index=sample_id, columns=[list(range(256))])
full_npy_data.to_csv(result_file+'emb.csv')


