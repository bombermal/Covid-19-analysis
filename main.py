#Imports
import pandas as pd

#My imports
import Codes.function as fc
import Codes.graph as gf
import Codes.read as rdF
#%% 1 - Load files

#Load trimmed aln
path = "F:/Ivan/OneDrive/Códigos ( Profissional )/IC/Covid/Covid-19-analysis/read/"
file = "gisaid_cov26032020_sequences_trimmed.aln.fasta"
loaded_files = rdF.read_aligned_files(path, file)
#Put Sequences in a dataframe
raw_df_aln = rdF.seq_to_df(loaded_files)
#Put sequences in a list, each row as a element
raw_lst_of_samples = fc.df_to_list(raw_df_aln)

#Result of processing od raw_df_aln, each row contains a (Position, Nucleotide)
final_data = pd.read_csv("Saved/Tidy_DF_final.csv")
poli_data = pd.read_csv("Saved/New_counting_df.csv")
#%% Divide o trabalho para poder ser executado em varios computadores
#Lista dos nomes dos arquivos que serão lidos
#divided_df_file_names = ["df_1.csv","df_2.csv","df_3.csv","df_4.csv"]
#DataFrame utilizado para a criação da tabela de contagem de nucleotídeos
#new_counting_df = pd.DataFrame(data=0, index=range(1,29355), columns=final_data.Nuc.unique() )


#%%

# for idx, row in tqdm(final_data.iterrows(), total=len(final_data)):
#     plus_one(new_counting_df, row)
    
# new_counting_df.to_csv("Saved/New_counting_df.csv")

#%% Ploting

conditions = [.5,1,2]
list_of_filtered_dfs = []

for ii in conditions:
  list_of_filtered_dfs.append(fc.filter_criteria(poli_data, 1302, ii))
  
gf.three_plots(conditions, list_of_filtered_dfs)