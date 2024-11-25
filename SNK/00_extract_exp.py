import os
import sys
import pandas as pd

path = sys.argv[1]
output = sys.argv[2]
tab_path_list = []
name_list = []

# 获取目标目录和文件
for root, dirs, files in os.walk(path, topdown=False):
    for name in files:
        name_tab = name.split('.')
        if name_tab[-1] == "tab":
            name_list.append(name_tab[0])
            tab_path_list.append(os.path.join(root, name))

sample_mapping = dict(zip(name_list,tab_path_list))

merge_pd = pd.DataFrame()
# 合并fpkm
num = 1
for tab_name, tab_path in sample_mapping.items():
    extracted_pd = pd.read_table(tab_path,sep='\t',header=0,index_col=0)
    if num == 1:
        merge_pd = pd.DataFrame(index = extracted_pd.index)
        merge_pd[tab_name] = extracted_pd['FPKM']
        print(tab_name) 
        num += 1
    else:
        merge_pd[tab_name] = extracted_pd['FPKM']
        print(tab_name)
merge_pd.sort_index(axis=1,inplace=True)
merge_pd.to_csv(output,sep='\t')
