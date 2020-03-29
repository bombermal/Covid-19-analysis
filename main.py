#Imports
import pandas as pd
import Codes.read as rdF
import seaborn as sns
from matplotlib import pyplot
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp

#%%

#Load files
path = "F:/Ivan/OneDrive/Códigos ( Profissional )/IC/Covid/Covid-19-analysis/read/"
file = "gisaid_cov26032020_sequences_trimmed.aln.fasta"
loaded_files = rdF.read_aligned_files(path, file)

#%%

#Save on DF
data_frame = rdF.seq_to_df(loaded_files)

#%%
def single_nuc_count(sample):
    aux = pd.DataFrame(columns=["Pos", "Nuc"])

    for idx, value in enumerate(sample.Seq):
        aux.loc[len(aux)] = [idx+1, value]  
    return aux

def df_to_list(df):
    temp = []
    for idx, row in df.iterrows():
        temp.append(row)
        
    return temp

list_of_samples = df_to_list(data_frame)

#%%
# Multiprocessing Imap
tidy_count_df2 = pd.DataFrame(columns=["Pos", "Nuc"])
poolSize = mp.cpu_count()
pool = ThreadPool(poolSize)

for ii in tqdm(pool.imap_unordered(single_nuc_count, list_of_samples), total=len(list_of_samples)):
    tidy_count_df2 = pd.concat([tidy_count_df2, ii])
pool.close()
pool.join()

tidy_count_df2.to_csv("Tidy_DF.csv", index=False)
#%%
a4_dims = (511.7, 8.27)
fig, ax = pyplot.subplots(figsize=a4_dims)
sns.countplot(ax=ax, x="Pos", hue="Nuc", data=tidy_count_df2)

#%%


