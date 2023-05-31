import pandas as pd
import numpy as np
from collections import Counter
import os


######################################################### gw11
os.chdir('/data3/yuhan/Project/lr/URD/yyc_result/')

# load raw transition and cell subtypes annotation
df = pd.read_csv('axialGW11.allTransition.txt', index_col=0, sep='\t')
df_anno = pd.read_csv('axialGW11.meta.txt', index_col=0, sep='\t')
dic_ct = dict(zip(df_anno.index, df_anno['cellType']))

df['from_ct'] = df['from'].map(dic_ct)
df['to_ct'] = df['to'].map(dic_ct)

# random select transitions between specific cell subtypes
df_retain = pd.concat([df[(df['from_ct']=='RG') & (df['to_ct']=='EOMES_IPC')].sample(n=500),
                        df[(df['from_ct']=='EOMES_IPC') & (df['to_ct']=='EN')].sample(n=700), 
                        df[(df['from_ct']=='RG') & (df['to_ct']=='ASCL1_IPC')].sample(n=500), 
                        df[(df['from_ct']=='ASCL1_IPC') & (df['to_ct']=='IN_Precursor')].sample(n=700), 
                        df[(df['from_ct']=='IN_Precursor') & (df['to_ct']=='IN')].sample(n=500),
                        df[(df['from_ct']=='IN') & (df['to_ct']=='IN')].sample(n=500),
                        df[(df['from_ct']=='EN') & (df['to_ct']=='EN')].sample(n=700)])


# random selecte 2 transitions for each cell
all_cells = set(df['from']) | set(df['to'])
random_index = []
a = 0
for i in all_cells:
    print(a, end='\r')
    a+=1
    index_tmp = df[(df['from'] == i) | (df['to'] == i)].index
    cur_tmp = np.random.choice(index_tmp, 2)
    random_index.extend(list(cur_tmp))


# remove some unexpected transitions
df_remove = pd.concat([df[(df['from_ct']=='IN') & (df['to_ct']=='EN')], df[(df['from_ct']=='EN') & (df['to_ct']=='IN')]])
remain_edge = set(df_retain.index) | set(random_index)
df_remain = df.reindex(remain_edge)
for i in df_remove.index:
    if i in df_remain.index:
        c1 = df_remain.loc[i, 'from']
        c2 = df_remain.loc[i, 'to']
        cur_c = Counter(list(df_remain['from']) + list(df_remain['to']))
        if cur_c[c1] > 1 and cur_c[c2] > 1:
            df_remain.drop(i, inplace=True)

df_remain.to_csv('axialGW11.UsedTransition.txt', sep='\t')



############################################################ GW20
# load raw transition and cell subtypes annotation
df = pd.read_csv('axialGW20.allTransition.txt', index_col=0, sep='\t')
df_anno = pd.read_csv('axialGW20.meta.txt', index_col=0, sep='\t')
dic_ct = dict(zip(df_anno.index, df_anno['cellType']))

df['from_ct'] = df['from'].map(dic_ct)
df['to_ct'] = df['to'].map(dic_ct)

# random select transitions between specific cell subtypes
df_retain = pd.concat([df[(df['from_ct']=='RG') & (df['to_ct']=='EOMES_IPC')].sample(n=400),
                        df[(df['from_ct']=='EOMES_IPC') & (df['to_ct']=='EN')].sample(n=1000), 
                        df[(df['from_ct']=='RG') & (df['to_ct']=='ASCL1_IPC')].sample(n=800), 
                        df[(df['from_ct']=='ASCL1_IPC') & (df['to_ct']=='IN_Precursor')].sample(n=1000), 
                        df[(df['from_ct']=='IN_Precursor') & (df['to_ct']=='IN')].sample(n=1000),
                        df[(df['from_ct']=='IN') & (df['to_ct']=='IN')].sample(n=1200),
                        df[(df['from_ct']=='EN') & (df['to_ct']=='EN')].sample(n=1000)])


# random selecte 2 transitions for each cell
all_cells = set(df['from']) | set(df['to'])
random_index = []
a = 0
for i in all_cells:
    print(a, end='\r')
    a+=1
    index_tmp = df[(df['from'] == i) | (df['to'] == i)].index
    cur_tmp = np.random.choice(index_tmp, 2)
    random_index.extend(list(cur_tmp))

# remove some unexpected transitions
df_remove = pd.concat([df[(df['from_ct']=='IN') & (df['to_ct']=='EN')], df[(df['from_ct']=='EN') & (df['to_ct']=='IN')]])
remain_edge = set(df_retain.index) | set(random_index)
df_remain = df.reindex(remain_edge)
for i in df_remove.index:
    if i in df_remain.index:
        c1 = df_remain.loc[i, 'from']
        c2 = df_remain.loc[i, 'to']
        cur_c = Counter(list(df_remain['from']) + list(df_remain['to']))
        if cur_c[c1] > 1 and cur_c[c2] > 1:
            df_remain.drop(i, inplace=True)
df_remain.to_csv('axialGW20.UsedTransition.txt', sep='\t')



################################################## merge
# merge gw11 and gw20 transitions
df1 = pd.read_csv('axialGW11.UsedTransition.txt', index_col=0, sep='\t')
df2 = pd.read_csv('axialGW20.UsedTransition.txt', index_col=0, sep='\t')

df = pd.concat((df1, df2))
df.index = range(df.shape[0])
df_remove = df.sample(n=15000)

for i in df_remove.index:
    c1 = df.loc[i, 'from']
    c2 = df.loc[i, 'to']
    cur_c = Counter(list(df['from']) + list(df['to']))
    if cur_c[c1] > 1 and cur_c[c2] > 1:
        df.drop(i, inplace=True)

df.to_csv('axialMerge.UsedTransition.txt', sep='\t')