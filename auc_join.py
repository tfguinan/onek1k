import os
import pandas as pd
import loompy as lp
from glob import glob
from pyscenic.rss import regulon_specificity_scores

in_path = '/data/menzies_projects/onek1k/share/TG/git_analysis/output_data/scenic/'
loom_list = glob(f'{in_path}*-scenic-output.loom')

# Init empty to concat with data
full_auc_mtx = pd.DataFrame()
cell_type_index = pd.DataFrame()

# We should ignore index
# , ignore_index=True, sort=False

for loom in loom_list:
    print('Accessing:', loom)
    lf = lp.connect(loom, mode='r+', validate=False)
    
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    cell_type = pd.DataFrame(lf.ca['predicted.celltype.l2'], index=lf.ca.CellID)

    cell_type_index = pd.concat([cell_type_index, cell_type], ignore_index=True, sort=False)
    full_auc_mtx = pd.concat([full_auc_mtx, auc_mtx], ignore_index=True, sort=False)

# Convert NaN to 0
full_auc_mtx = full_auc_mtx.fillna(0)

cell_type_index.columns = ['predicted.celltype.l2']
print(cell_type_index)

# HUGE! ~ 14GB
# print('Writing to csv')
# full_auc_mtx.to_csv('/data/menzies_projects/onek1k/share/TG/git_analysis/full_auc_mtx.csv')

print('Calculating RSS')
rss_cell_type = regulon_specificity_scores(full_auc_mtx, cell_type_index['predicted.celltype.l2'] )

print('Writing to csv')
rss_cell_type.to_csv('/data/menzies_projects/onek1k/share/TG/git_analysis/cell_type_rss.csv')